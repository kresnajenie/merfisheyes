"use client";

import {
  Navbar as HeroUINavbar,
  NavbarContent,
  NavbarMenu,
  NavbarMenuToggle,
  NavbarBrand,
  NavbarItem,
  NavbarMenuItem,
} from "@heroui/navbar";
import { Button } from "@heroui/button";
import { Kbd } from "@heroui/kbd";
import { Link } from "@heroui/link";
import { Input } from "@heroui/input";
import { link as linkStyles } from "@heroui/theme";
import NextLink from "next/link";
import clsx from "clsx";
import { usePathname } from "next/navigation";
import { useEffect, useState } from "react";

import { siteConfig } from "@/config/site";
import {
  TwitterIcon,
  GithubIcon,
  DiscordIcon,
  HeartFilledIcon,
  SearchIcon,
} from "@/components/icons";
import { useDatasetStore } from "@/lib/stores/datasetStore";
import { useSingleMoleculeStore } from "@/lib/stores/singleMoleculeStore";

interface NavbarProps {
  onUploadClick?: () => void;
}

export const Navbar = ({ onUploadClick }: NavbarProps) => {
  const pathname = usePathname();
  const cellDatasetId = useDatasetStore((state) => state.currentDatasetId);
  const smDatasetId = useSingleMoleculeStore(
    (state) => state.currentDatasetId,
  );
  const [viewerHref, setViewerHref] = useState(() => {
    if (typeof window !== "undefined") {
      const storedMode = window.localStorage.getItem("lastDatasetMode");

      if (storedMode === "sm") {
        return "/sm-viewer";
      }
    }

    return "/viewer";
  });

  useEffect(() => {
    if (typeof window === "undefined") return;

    const storedMode = window.localStorage.getItem("lastDatasetMode");

    if (storedMode === "sm") {
      setViewerHref("/sm-viewer");

      return;
    }

    if (storedMode === "cell") {
      setViewerHref("/viewer");

      return;
    }

    if (smDatasetId && !cellDatasetId) {
      setViewerHref("/sm-viewer");
    } else {
      setViewerHref("/viewer");
    }
  }, [cellDatasetId, smDatasetId]);

  const searchInput = (
    <Input
      aria-label="Search"
      classNames={{
        inputWrapper: "bg-default-100",
        input: "text-sm",
      }}
      endContent={
        <Kbd className="hidden lg:inline-block" keys={["command"]}>
          K
        </Kbd>
      }
      labelPlacement="outside"
      placeholder="Search..."
      startContent={
        <SearchIcon className="text-base text-default-400 pointer-events-none flex-shrink-0" />
      }
      type="search"
    />
  );

  return (
    <HeroUINavbar
      className="w-3/4 mx-auto my-4 rounded-full glass"
      maxWidth="xl"
      position="sticky"
    >
      <NavbarContent className="basis-1/5 sm:basis-full" justify="start">
        <NavbarBrand as="li" className="gap-3 max-w-fit">
          <NextLink className="flex justify-start items-center gap-1" href="/">
            {/* <Logo /> */}
            <p className="font-bold text-inherit">MERFISHEYES</p>
          </NextLink>
        </NavbarBrand>
        <ul className="hidden lg:flex gap-4 justify-start ml-2">
          {siteConfig.navItems.map((item) => {
            const href = item.href === "/viewer" ? viewerHref : item.href;

            return (
              <NavbarItem key={item.href}>
                <NextLink
                  className={clsx(
                    linkStyles({ color: "foreground" }),
                    "data-[active=true]:text-primary data-[active=true]:font-medium",
                  )}
                  color="foreground"
                  href={href}
                >
                  {item.label}
                </NextLink>
              </NavbarItem>
            );
          })}
        </ul>
      </NavbarContent>

      <NavbarContent
        className="hidden sm:flex basis-1/5 sm:basis-full"
        justify="end"
      >
        {/* Show Upload & Save button on /viewer when dataset is loaded */}
        {pathname === "/viewer" && onUploadClick && (
          <NavbarItem>
            <Button color="primary" onPress={onUploadClick}>
              Upload & Save
            </Button>
          </NavbarItem>
        )}

        {/* Show Upload & Save button on /sm-viewer when dataset is loaded */}
        {pathname === "/sm-viewer" && onUploadClick && (
          <NavbarItem>
            <Button color="primary" onPress={onUploadClick}>
              Upload & Save
            </Button>
          </NavbarItem>
        )}

        {/* Show Back to Viewer button on /viewer/[id] */}
        {pathname?.startsWith("/viewer/") && pathname !== "/viewer" && (
          <NavbarItem>
            <Button as={NextLink} color="primary" href="/viewer" variant="flat">
              Back to Viewer
            </Button>
          </NavbarItem>
        )}

        <NavbarItem className="hidden sm:flex gap-2">
          <Link isExternal aria-label="Twitter" href={siteConfig.links.twitter}>
            <TwitterIcon className="text-default-500" />
          </Link>
          <Link isExternal aria-label="Discord" href={siteConfig.links.discord}>
            <DiscordIcon className="text-default-500" />
          </Link>
          <Link isExternal aria-label="Github" href={siteConfig.links.github}>
            <GithubIcon className="text-default-500" />
          </Link>
        </NavbarItem>
        {/* <NavbarItem className="hidden lg:flex">{searchInput}</NavbarItem> */}
        <NavbarItem className="hidden md:flex">
          <Button
            isExternal
            as={Link}
            className="text-sm font-normal text-default-600 bg-default-100"
            href={siteConfig.links.sponsor}
            startContent={<HeartFilledIcon className="text-danger" />}
            variant="flat"
          >
            Sponsor
          </Button>
        </NavbarItem>
      </NavbarContent>

      <NavbarContent className="sm:hidden basis-1 pl-4" justify="end">
        <Link isExternal aria-label="Github" href={siteConfig.links.github}>
          <GithubIcon className="text-default-500" />
        </Link>
        <NavbarMenuToggle />
      </NavbarContent>

      <NavbarMenu>
        {searchInput}
        <div className="mx-4 mt-2 flex flex-col gap-2">
          {siteConfig.navMenuItems.map((item, index) => {
            const href = item.href === "/viewer" ? viewerHref : item.href;

            return (
              <NavbarMenuItem key={`${item}-${index}`}>
                <Link color="foreground" href={href} size="lg">
                  {item.label}
                </Link>
              </NavbarMenuItem>
            );
          })}
        </div>
      </NavbarMenu>
    </HeroUINavbar>
  );
};
