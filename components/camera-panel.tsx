"use client";

import { Button } from "@heroui/button";
import { Slider } from "@heroui/react";
import { useRef, useEffect } from "react";

import { glassButton } from "@/components/primitives";

interface CameraPanelProps {
  onClose: () => void;
  controlsRef?: React.RefObject<HTMLDivElement>;
  sceneRotation: number;
  setSceneRotation: (degrees: number) => void;
  flipX: boolean;
  setFlipX: (flip: boolean) => void;
  flipY: boolean;
  setFlipY: (flip: boolean) => void;
  viewMode?: string;
  setViewMode?: (mode: any) => void;
  is3DDataset?: boolean;
}

export function CameraPanel({
  onClose,
  controlsRef,
  sceneRotation,
  setSceneRotation,
  flipX,
  setFlipX,
  flipY,
  setFlipY,
  viewMode,
  setViewMode,
  is3DDataset = false,
}: CameraPanelProps) {
  const panelRef = useRef<HTMLDivElement>(null);

  // Handle click outside to close panel
  useEffect(() => {
    const handleClickOutside = (event: MouseEvent) => {
      const target = event.target as Node;

      if (panelRef.current && panelRef.current.contains(target)) return;
      if (controlsRef?.current && controlsRef.current.contains(target)) return;

      onClose();
    };

    document.addEventListener("mousedown", handleClickOutside);
    return () => document.removeEventListener("mousedown", handleClickOutside);
  }, [onClose, controlsRef]);

  return (
    <div
      ref={panelRef}
      className={`absolute top-0 left-16 z-50 w-[240px] border-2 border-white/20 rounded-3xl shadow-lg ${glassButton()}`}
    >
      <div className="p-4 space-y-4">
        {/* Header */}
        <div className="flex justify-between items-center">
          <span className="text-sm font-medium">Camera</span>
          <Button
            className="h-6 w-6 min-w-0"
            size="sm"
            variant="light"
            onPress={onClose}
          >
            <svg className="w-3 h-3" fill="none" stroke="currentColor" strokeWidth={2} viewBox="0 0 24 24">
              <path d="M6 18L18 6M6 6l12 12" strokeLinecap="round" />
            </svg>
          </Button>
        </div>

        {/* Rotation */}
        <div>
          <div className="flex justify-between text-xs mb-1">
            <span className="text-default-500">Rotation</span>
            <span className="text-default-400">{sceneRotation}°</span>
          </div>
          <Slider
            aria-label="Rotation"
            className="w-full"
            maxValue={360}
            minValue={0}
            size="sm"
            step={1}
            value={sceneRotation}
            onChange={(v) => {
              let val = v as number;
              if ((window as any).__shiftHeld) {
                val = Math.round(val / 45) * 45;
              }
              setSceneRotation(val);
            }}
          />
          <p className="text-[10px] text-default-400 mt-1">Hold Shift for 45° snapping</p>
        </div>

        {/* Flip */}
        <div>
          <span className="text-xs text-default-500 mb-2 block">Flip</span>
          <div className="flex gap-2">
            <Button
              className="flex-1 text-xs"
              color={flipX ? "primary" : "default"}
              size="sm"
              variant={flipX ? "solid" : "bordered"}
              onPress={() => setFlipX(!flipX)}
            >
              ⇔ Horizontal
            </Button>
            <Button
              className="flex-1 text-xs"
              color={flipY ? "primary" : "default"}
              size="sm"
              variant={flipY ? "solid" : "bordered"}
              onPress={() => setFlipY(!flipY)}
            >
              ⇕ Vertical
            </Button>
          </div>
        </div>

        {/* 2D/3D Toggle */}
        {is3DDataset && setViewMode && viewMode && (
          <div>
            <span className="text-xs text-default-500 mb-2 block">View Mode</span>
            <div className="flex gap-2">
              <Button
                className="flex-1 text-xs"
                color={viewMode === "2D" ? "primary" : "default"}
                size="sm"
                variant={viewMode === "2D" ? "solid" : "bordered"}
                onPress={() => setViewMode("2D")}
              >
                2D
              </Button>
              <Button
                className="flex-1 text-xs"
                color={viewMode === "3D" ? "primary" : "default"}
                size="sm"
                variant={viewMode === "3D" ? "solid" : "bordered"}
                onPress={() => setViewMode("3D")}
              >
                3D
              </Button>
            </div>
          </div>
        )}

        {/* Reset */}
        <Button
          className="w-full"
          color="danger"
          size="sm"
          variant="flat"
          onPress={() => {
            setSceneRotation(0);
            setFlipX(false);
            setFlipY(false);
          }}
        >
          Reset
        </Button>
      </div>
    </div>
  );
}
