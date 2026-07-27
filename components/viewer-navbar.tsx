"use client";

import { Button } from "@heroui/button";
import { Link } from "@heroui/link";
import { link as linkStyles } from "@heroui/theme";
import NextLink from "next/link";
import clsx from "clsx";
import { usePathname } from "next/navigation";

import { siteConfig } from "@/config/site";
import {
  TwitterIcon,
  GithubIcon,
  DiscordIcon,
  HeartFilledIcon,
} from "@/components/icons";
import { glassButton } from "@/components/primitives";
import { UserMenu } from "@/components/auth/user-menu";

interface ViewerNavbarProps {
  onUploadClick?: () => void;
}

/**
 * Compact top-left brand pill for the viewers. Collapsed it shows only
 * "MERFISHEYES" (which links home, same as the logo elsewhere); on hover it
 * expands horizontally to reveal the full nav — links, socials, Sponsor, and
 * the Upload & Save action on from-local viewers.
 */
export const ViewerNavbar = ({ onUploadClick }: ViewerNavbarProps) => {
  const pathname = usePathname();
  const showUpload =
    (pathname === "/viewer/from-local" ||
      pathname === "/sm-viewer/from-local") &&
    !!onUploadClick;

  return (
    <div className="absolute top-4 left-4 z-[var(--z-chrome)] flex items-center gap-2">
      <div
        className={clsx(
          "group flex items-center h-12 pl-4 pr-4 rounded-full overflow-hidden",
          glassButton(),
        )}
      >
        {/* Brand — always visible, links home */}
        <NextLink
          className="font-bold text-white whitespace-nowrap shrink-0"
          href="/"
        >
          MERFISHEYES
        </NextLink>

        {/* Everything else — collapses to zero width until the pill is hovered */}
        <div className="grid grid-cols-[0fr] group-hover:grid-cols-[1fr] transition-[grid-template-columns] duration-300 ease-out">
          <div className="overflow-hidden">
            <div className="flex items-center gap-4 pl-4 whitespace-nowrap opacity-0 group-hover:opacity-100 transition-opacity duration-300">
              {siteConfig.navItems.map((item) => (
                <NextLink
                  key={item.href}
                  className={clsx(linkStyles({ color: "foreground" }))}
                  color="foreground"
                  href={item.href}
                >
                  {item.label}
                </NextLink>
              ))}

              <div className="flex items-center gap-2">
                <Link
                  isExternal
                  aria-label="Twitter"
                  href={siteConfig.links.twitter}
                >
                  <TwitterIcon className="text-default-500" />
                </Link>
                <Link
                  isExternal
                  aria-label="Discord"
                  href={siteConfig.links.discord}
                >
                  <DiscordIcon className="text-default-500" />
                </Link>
                <Link
                  isExternal
                  aria-label="Github"
                  href={siteConfig.links.github}
                >
                  <GithubIcon className="text-default-500" />
                </Link>
              </div>

              <Button
                isExternal
                as={Link}
                className="text-sm font-normal text-default-600 bg-default-100 shrink-0"
                href={siteConfig.links.sponsor}
                size="sm"
                startContent={<HeartFilledIcon className="text-danger" />}
                variant="flat"
              >
                Sponsor
              </Button>

              <UserMenu />
            </div>
          </div>
        </div>
      </div>

      {/* Upload & Save — a primary action for locally loaded datasets, kept
          always visible (not behind the hover) so it's a one-click save. */}
      {showUpload && (
        <Button color="primary" onPress={onUploadClick}>
          Upload & Save
        </Button>
      )}
    </div>
  );
};
