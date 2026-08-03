"use client";

import { Button } from "@heroui/button";
import {
  Modal,
  ModalBody,
  ModalContent,
  ModalFooter,
  ModalHeader,
  Textarea,
  useDisclosure,
} from "@heroui/react";
import { useEffect, useRef, useState } from "react";
import { toast } from "react-toastify";

import { usePanelVisualizationStore } from "@/lib/hooks/usePanelStores";
import { useTourStore } from "@/lib/stores/tourStore";
import { isValidTourConfig, type TourConfig } from "@/lib/tour/types";

const EMPTY_TOUR_TEMPLATE = `{
  "defaults": {
    "flyInMs": 1000,
    "orbitMs": 1000,
    "orbitRevolutions": 1,
    "zoomOutMs": 1000,
    "zoomRadius": 300,
    "zoomOutRadius": 1500,
    "filterRadius": 30
  },
  "alwaysOnGenes": [
    { "name": "GENE_BG", "color": "#ff0000" }
  ],
  "stops": [
    { "name": "Stop 1", "x": 0, "y": 0, "z": 0, "genes": ["GENE_A"] }
  ]
}
`;

export function TourSection() {
  const fileInputRef = useRef<HTMLInputElement>(null);
  const config = useTourStore((s) => s.config);
  const configName = useTourStore((s) => s.configName);
  const isPlaying = useTourStore((s) => s.isPlaying);
  const currentStopIndex = useTourStore((s) => s.currentStopIndex);
  const currentPhase = useTourStore((s) => s.currentPhase);
  const setConfig = useTourStore((s) => s.setConfig);
  const start = useTourStore((s) => s.start);
  const stop = useTourStore((s) => s.stop);
  const viewMode = usePanelVisualizationStore((s) => s.viewMode);

  const editor = useDisclosure();

  const handleFile = (file: File) => {
    const reader = new FileReader();

    reader.onload = () => {
      try {
        const parsed = JSON.parse(reader.result as string);

        if (!isValidTourConfig(parsed)) {
          toast.error(
            "Invalid tour file: expected { stops: [{ name, x, y, z, genes }] }",
          );

          return;
        }
        setConfig(parsed, file.name);
        toast.success(`Loaded tour: ${parsed.stops.length} stops`);
      } catch (err) {
        toast.error(`Failed to parse tour JSON: ${(err as Error).message}`);
      }
    };
    reader.onerror = () => toast.error("Failed to read tour file");
    reader.readAsText(file);
  };

  const handlePlay = () => {
    if (viewMode !== "3D") {
      toast.warning("Switch to 3D view before playing the tour");

      return;
    }
    start();
  };

  const currentStopName =
    config && isPlaying ? (config.stops[currentStopIndex]?.name ?? "") : "";

  const handleDownload = () => {
    if (!config) return;
    const blob = new Blob([JSON.stringify(config, null, 2)], {
      type: "application/json",
    });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");

    a.href = url;
    a.download = configName ?? "tour.json";
    a.click();
    URL.revokeObjectURL(url);
  };

  return (
    <div className="space-y-2">
      <span className="text-xs font-medium text-primary">Tour</span>

      <input
        ref={fileInputRef}
        accept="application/json,.json"
        className="hidden"
        type="file"
        onChange={(e) => {
          const f = e.target.files?.[0];

          if (f) handleFile(f);
          e.target.value = "";
        }}
      />

      <div className="flex gap-1">
        <Button
          className="flex-1"
          size="sm"
          variant="flat"
          onPress={() => fileInputRef.current?.click()}
        >
          {config ? "Replace" : "Load"}
        </Button>
        <Button
          className="flex-1"
          size="sm"
          variant="flat"
          onPress={editor.onOpen}
        >
          Edit
        </Button>
        <Button
          className="flex-1"
          isDisabled={!config}
          size="sm"
          variant="flat"
          onPress={handleDownload}
        >
          Save
        </Button>
      </div>

      {config && (
        <div className="text-xs text-default-500 space-y-0.5">
          <div className="truncate" title={configName ?? undefined}>
            {configName ?? "tour.json"}
          </div>
          <div>
            {config.stops.length} stop{config.stops.length === 1 ? "" : "s"}
          </div>
          {isPlaying && (
            <div className="text-default-400">
              {currentStopIndex + 1}/{config.stops.length} · {currentPhase} ·{" "}
              <span className="text-default-600">{currentStopName}</span>
            </div>
          )}
        </div>
      )}

      <Button
        className="w-full"
        color={isPlaying ? "danger" : "primary"}
        isDisabled={!config}
        size="sm"
        variant="flat"
        onPress={() => (isPlaying ? stop() : handlePlay())}
      >
        {isPlaying ? "Stop tour" : "Play tour"}
      </Button>

      <TourJsonEditor
        initialName={configName}
        initialText={
          config ? JSON.stringify(config, null, 2) : EMPTY_TOUR_TEMPLATE
        }
        isOpen={editor.isOpen}
        onApply={(parsed, name) => {
          setConfig(parsed, name ?? configName ?? "tour.json");
          toast.success(`Saved tour: ${parsed.stops.length} stops`);
          editor.onClose();
        }}
        onClose={editor.onClose}
      />
    </div>
  );
}

interface TourJsonEditorProps {
  isOpen: boolean;
  onClose: () => void;
  initialText: string;
  initialName: string | null;
  onApply: (parsed: TourConfig, name: string | null) => void;
}

function TourJsonEditor({
  isOpen,
  onClose,
  initialText,
  initialName,
  onApply,
}: TourJsonEditorProps) {
  const [text, setText] = useState(initialText);
  const [name, setName] = useState(initialName ?? "tour.json");
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    if (isOpen) {
      setText(initialText);
      setName(initialName ?? "tour.json");
      setError(null);
    }
  }, [isOpen, initialText, initialName]);

  const handleApply = () => {
    let parsed: unknown;

    try {
      parsed = JSON.parse(text);
    } catch (err) {
      setError(`JSON parse error: ${(err as Error).message}`);

      return;
    }
    if (!isValidTourConfig(parsed)) {
      setError(
        "Schema invalid: expected { stops: [{ name, x, y, z, genes: [...] }] }",
      );

      return;
    }
    setError(null);
    onApply(parsed, name);
  };

  return (
    <Modal isOpen={isOpen} scrollBehavior="inside" size="3xl" onClose={onClose}>
      <ModalContent>
        <ModalHeader className="flex flex-col gap-1">
          <span>Edit tour JSON</span>
          <span className="text-xs font-normal text-default-500">
            Changes apply immediately when you click Apply. Use Download in the
            panel to save the result to disk.
          </span>
        </ModalHeader>
        <ModalBody>
          <input
            className="w-full bg-default-100 rounded px-2 py-1 text-sm"
            placeholder="filename.json"
            type="text"
            value={name}
            onChange={(e) => setName(e.target.value)}
          />
          <Textarea
            classNames={{
              input: "font-mono text-xs leading-snug",
            }}
            minRows={20}
            value={text}
            onValueChange={setText}
          />
          {error && (
            <div className="text-xs text-danger whitespace-pre-wrap">
              {error}
            </div>
          )}
        </ModalBody>
        <ModalFooter>
          <Button size="sm" variant="light" onPress={onClose}>
            Cancel
          </Button>
          <Button color="primary" size="sm" onPress={handleApply}>
            Apply
          </Button>
        </ModalFooter>
      </ModalContent>
    </Modal>
  );
}
