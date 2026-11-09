import { type ChildProcessWithoutNullStreams, spawn } from "node:child_process";
import path from "node:path";

/**
 * Drive the built merfisheyes-mcp package over stdio as an MCP client would
 * — newline-delimited JSON-RPC (MCP stdio transport) — so a test sees exactly
 * what an agent sees. No MCP SDK needed on the site's side.
 */

export interface ToolResult {
  isError?: boolean;
  content: { type: string; text: string }[];
  structuredContent?: Record<string, unknown>;
}

export class McpStdio {
  private next = 1;
  private waiting = new Map<number, (msg: { result?: unknown; error?: { message: string } }) => void>();
  private buffer = "";
  stderr = "";

  private constructor(private proc: ChildProcessWithoutNullStreams) {
    proc.stdout.setEncoding("utf8").on("data", (chunk: string) => {
      this.buffer += chunk;
      let nl: number;

      while ((nl = this.buffer.indexOf("\n")) >= 0) {
        const line = this.buffer.slice(0, nl).trim();

        this.buffer = this.buffer.slice(nl + 1);
        if (!line) continue;
        const msg = JSON.parse(line);

        this.waiting.get(msg.id)?.(msg);
        this.waiting.delete(msg.id);
      }
    });
    proc.stderr.setEncoding("utf8").on("data", (c: string) => (this.stderr += c));
  }

  /** Start the package (packages/merfisheyes-mcp/dist) and initialize. */
  static async start(env: Record<string, string>): Promise<McpStdio> {
    const entry = path.join(__dirname, "../../packages/merfisheyes-mcp/dist/index.js");
    const proc = spawn(process.execPath, [entry], { env: { PATH: process.env.PATH ?? "", HOME: process.env.HOME ?? "", ...env } });
    const client = new McpStdio(proc);

    await client.request("initialize", {
      protocolVersion: "2025-06-18",
      capabilities: {},
      clientInfo: { name: "merfisheyes-e2e", version: "1.0.0" },
    });
    client.notify("notifications/initialized");

    return client;
  }

  request<T = unknown>(method: string, params: unknown = {}): Promise<T> {
    const id = this.next++;

    return new Promise<T>((resolve, reject) => {
      const timer = setTimeout(() => reject(new Error(`${method} timed out; stderr:\n${this.stderr}`)), 30_000);

      this.waiting.set(id, (msg) => {
        clearTimeout(timer);
        if (msg.error) reject(new Error(`${method}: ${msg.error.message}`));
        else resolve(msg.result as T);
      });
      this.proc.stdin.write(`${JSON.stringify({ jsonrpc: "2.0", id, method, params })}\n`);
    });
  }

  notify(method: string, params: unknown = {}) {
    this.proc.stdin.write(`${JSON.stringify({ jsonrpc: "2.0", method, params })}\n`);
  }

  call(name: string, args: Record<string, unknown> = {}): Promise<ToolResult> {
    return this.request<ToolResult>("tools/call", { name, arguments: args });
  }

  close() {
    this.proc.kill();
  }
}
