import "@/styles/globals.css";
import { Metadata, Viewport } from "next";
import { Link } from "@heroui/link";
import clsx from "clsx";

import { Providers } from "./providers";

import { siteConfig } from "@/config/site";
import { fontSans } from "@/config/fonts";
import { Navbar } from "@/components/navbar";

export const metadata: Metadata = {
  title: {
    default: siteConfig.name,
    template: `%s - ${siteConfig.name}`,
  },
  description: siteConfig.description,
  icons: {
    icon: "/favicon.ico",
  },
};

export const viewport: Viewport = {
  themeColor: [
    { media: "(prefers-color-scheme: light)", color: "white" },
    { media: "(prefers-color-scheme: dark)", color: "black" },
  ],
};

export default function RootLayout({
  children,
}: {
  children: React.ReactNode;
}) {
  return (
    <html suppressHydrationWarning lang="en">
      <head />
      <body
        className={clsx(
          "min-h-screen text-foreground bg-background font-sans antialiased",
          fontSans.variable
        )}
      >
        <Providers themeProps={{ attribute: "class", defaultTheme: "dark" }}>
          <div className="relative flex flex-col h-screen">
            <Navbar />
            <main className="container mx-auto max-w-9xl pt-16 px-6 flex-grow">
              {children}
            </main>
            <footer className="w-full flex items-center justify-center py-3">
              <Link
                isExternal
                className="flex items-center gap-1 text-current"
                href="https://www.linkedin.com/in/ignatius-jenie-1023521b3/"
                title="linkedin.com kresnajenie"
              >
                <span className="text-default-600">Made by</span>
                <p className="text-primary">Kresna Jenie</p>
              </Link>
              <Link
                isExternal
                className="flex items-center gap-1 text-current"
                href="https://b.bintulab.com"
                title="b.bintulab.com"
              >
                <span className="text-default-600"></span>
                <span className="text-default-600">at</span>
                <p className="text-primary">Bintu Lab</p>
              </Link>
            </footer>
          </div>
        </Providers>
      </body>
    </html>
  );
}
