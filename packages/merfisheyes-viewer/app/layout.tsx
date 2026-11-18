import "@/styles/globals.css";
import clsx from "clsx";

import { Providers } from "@/app/providers";
import { AgentSession } from "@/components/agent-session";
import { fontSans } from "@/config/fonts";

/** The standalone viewer's frame: the site's styles and providers, no
 * navbar, footer, analytics or account. */
export default function StandaloneLayout({
  children,
}: {
  children: React.ReactNode;
}) {
  return (
    <html suppressHydrationWarning lang="en">
      <head>
        <title>MERFISHEYES (local)</title>
      </head>
      <body
        className={clsx(
          "min-h-screen text-foreground bg-background font-sans antialiased",
          fontSans.variable,
        )}
      >
        <Providers
          noAccount
          themeProps={{ attribute: "class", defaultTheme: "dark" }}
        >
          <main className="h-screen w-screen">{children}</main>
          <AgentSession />
        </Providers>
      </body>
    </html>
  );
}
