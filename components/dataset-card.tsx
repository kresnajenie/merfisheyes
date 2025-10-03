"use client";

import { Card, CardBody, Image } from "@heroui/react";
import { useRouter } from "next/navigation";

export interface Dataset {
  id: string;
  title: string;
  description: string;
  image: string;
  type: string;
  link?: string;
}

interface DatasetCardProps {
  dataset: Dataset;
}

export function DatasetCard({ dataset }: DatasetCardProps) {
  const handlePress = () => {
    if (dataset.link) {
      window.open(dataset.link, "_blank", "noopener,noreferrer");
    }
  };

  return (
    <Card
      isBlurred
      isPressable
      onPress={handlePress}
      className="border-none bg-background/60 dark:bg-default-100/50 hover:bg-default-200/50 hover:scale-105 transition-all duration-200 cursor-pointer"
      shadow="sm"
    >
      <CardBody className="p-4">
        <div className="grid grid-cols-12 gap-4 items-center">
          {/* Image on the left */}
          <div className="relative col-span-4">
            <Image
              alt={dataset.title}
              className="object-cover"
              height={80}
              shadow="sm"
              src={dataset.image}
              width="100%"
            />
          </div>

          {/* Content on the right */}
          <div className="flex flex-col col-span-8">
            <div className="flex flex-col gap-1">
              <h3 className="font-semibold text-foreground text-sm">
                {dataset.title}
              </h3>
              <p className="text-xs text-foreground/60">{dataset.type}</p>
              <p className="text-xs text-foreground/80 mt-1">
                {dataset.description}
              </p>
            </div>
          </div>
        </div>
      </CardBody>
    </Card>
  );
}
