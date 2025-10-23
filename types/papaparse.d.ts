declare module "papaparse" {
  export type ParseInput = string | File | Blob;

  export interface ParseMeta {
    fields?: (string | undefined)[];
  }

  export interface ParseError {
    type: string;
    code: string;
    message: string;
    row?: number;
  }

  export interface ParseResult<T> {
    data: T[];
    errors: ParseError[];
    meta: ParseMeta;
  }

  export interface ParseParser {
    abort(): void;
    pause(): void;
    resume(): void;
  }

  export interface ParseConfig<T> {
    header?: boolean;
    skipEmptyLines?: boolean | "greedy";
    worker?: boolean;
    dynamicTyping?: boolean;
    transformHeader?: (header: string | undefined) => string;
    chunk?: (results: ParseResult<T>, parser: ParseParser) => void;
    complete?: (results: ParseResult<T>, file?: File) => void;
    error?: (error: unknown, file?: File) => void;
  }

  export interface PapaParser {
    // Synchronous parsing (string input without callbacks)
    parse<T>(input: string, config?: ParseConfig<T>): ParseResult<T>;
    // Asynchronous parsing (File/Blob input or with callbacks)
    parse<T>(input: File | Blob, config: ParseConfig<T>): void;
  }

  const Papa: PapaParser;
  export {
    ParseResult,
    ParseMeta,
    ParseError,
    ParseConfig,
    ParseParser,
    ParseInput,
  };
  export default Papa;
}
