"use client";

import { useEffect, useState } from "react";
import { Button } from "@heroui/button";
import NextLink from "next/link";

/**
 * Owner's Edit button on public Explore detail pages. The pages render
 * without auth() (for cacheability), so this resolves the viewer's edit
 * target client-side and renders nothing for everyone else.
 */
export function ExploreEditButton({ catalogId }: { catalogId: string }) {
  const [href, setHref] = useState<string | null>(null);

  useEffect(() => {
    let active = true;

    fetch(`/api/explore/${catalogId}/edit-target`)
      .then((res) => (res.ok ? res.json() : { editHref: null }))
      .then((j) => {
        if (active) setHref(j.editHref ?? null);
      })
      .catch(() => {});

    return () => {
      active = false;
    };
  }, [catalogId]);

  if (!href) return null;

  return (
    <Button as={NextLink} href={href} size="sm" variant="flat">
      Edit
    </Button>
  );
}
