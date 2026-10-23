"use client";

import type { ThemeProviderProps } from "next-themes";

import * as React from "react";
import { HeroUIProvider } from "@heroui/system";
import { SessionProvider } from "next-auth/react";
import { useRouter } from "next/navigation";
import { ThemeProvider as NextThemesProvider } from "next-themes";
import { ToastContainer, toast } from "react-toastify";
import "react-toastify/dist/ReactToastify.css";

import { HideUiManager } from "@/components/hide-ui-manager";
import { keepErrorsOpen } from "@/lib/utils/keep-errors-open";

// An error stays on screen until it is closed (lib/utils/keep-errors-open.ts)
keepErrorsOpen(toast);

export interface ProvidersProps {
  children: React.ReactNode;
  themeProps?: ThemeProviderProps;
}

declare module "@react-types/shared" {
  interface RouterConfig {
    routerOptions: NonNullable<
      Parameters<ReturnType<typeof useRouter>["push"]>[1]
    >;
  }
}

export function Providers({ children, themeProps }: ProvidersProps) {
  const router = useRouter();

  return (
    <SessionProvider>
      <HeroUIProvider navigate={router.push}>
        <NextThemesProvider {...themeProps}>
          {children}
          <HideUiManager />
          {/* A toast is not "outside" an open dialog: closing one (an error
              stays until it is) must not dismiss the dialog behind it, and
              react-aria keeps it readable to screen readers while a dialog
              is open */}
          <div data-react-aria-top-layer="true">
            <ToastContainer
              closeOnClick
              draggable
              pauseOnFocusLoss
              pauseOnHover
              autoClose={3000}
              hideProgressBar={false}
              newestOnTop={true}
              position="bottom-right"
              rtl={false}
              theme="dark"
            />
          </div>
        </NextThemesProvider>
      </HeroUIProvider>
    </SessionProvider>
  );
}
