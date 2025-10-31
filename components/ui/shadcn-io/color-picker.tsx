'use client';

import Color from 'color';
import { Pipette } from 'lucide-react';
import {
  type ComponentProps,
  createContext,
  type HTMLAttributes,
  useCallback,
  useContext,
  useEffect,
  useMemo,
  useRef,
  useState,
} from 'react';
import { Button } from '@heroui/button';
import { Input } from '@heroui/input';
import { Slider } from '@heroui/slider';
import {
  Select,
  SelectItem,
} from '@heroui/select';
import { cn } from '@/lib/utils';

interface ColorPickerContextValue {
  hue: number;
  saturation: number;
  lightness: number;
  alpha: number;
  mode: string;
  setHue: (hue: number) => void;
  setSaturation: (saturation: number) => void;
  setLightness: (lightness: number) => void;
  setAlpha: (alpha: number) => void;
  setMode: (mode: string) => void;
}

const ColorPickerContext = createContext<ColorPickerContextValue | undefined>(
  undefined
);

export const useColorPicker = () => {
  const context = useContext(ColorPickerContext);
  if (!context) {
    throw new Error('useColorPicker must be used within a ColorPickerProvider');
  }
  return context;
};

export type ColorPickerProps = Omit<HTMLAttributes<HTMLDivElement>, 'onChange'> & {
  value?: Parameters<typeof Color>[0];
  defaultValue?: Parameters<typeof Color>[0];
  onChange?: (value: string) => void;
};

export const ColorPicker = ({
  value,
  defaultValue = '#000000',
  onChange,
  className,
  children,
  ...props
}: ColorPickerProps) => {
  const initialColor = value ? Color(value) : Color(defaultValue);

  const [hue, setHue] = useState(initialColor.hue());
  const [saturation, setSaturation] = useState(initialColor.saturationl());
  const [lightness, setLightness] = useState(initialColor.lightness());
  const [alpha, setAlpha] = useState(1.0); // Always 100% opacity
  const [mode, setMode] = useState('hsl');
  const onChangeRef = useRef(onChange);

  // Keep ref updated
  useEffect(() => {
    onChangeRef.current = onChange;
  }, [onChange]);

  // Call onChange when internal state changes (not on mount)
  const isFirstRender = useRef(true);
  useEffect(() => {
    if (isFirstRender.current) {
      isFirstRender.current = false;
      return;
    }

    const color = Color.hsl(hue, saturation, lightness);
    // Use hex() instead of hexa() to exclude alpha channel
    onChangeRef.current?.(color.hex());
  }, [hue, saturation, lightness]);

  const contextValue = useMemo<ColorPickerContextValue>(
    () => ({
      hue,
      saturation,
      lightness,
      alpha,
      mode,
      setHue,
      setSaturation,
      setLightness,
      setAlpha,
      setMode,
    }),
    [hue, saturation, lightness, alpha, mode]
  );

  return (
    <ColorPickerContext.Provider value={contextValue}>
      <div className={cn('space-y-3', className)} {...props}>
        {children}
      </div>
    </ColorPickerContext.Provider>
  );
};

export const ColorPickerSelection = () => {
  const { hue, saturation, lightness } = useColorPicker();
  const [pos, setPos] = useState({ x: saturation, y: 100 - lightness });
  const canvasRef = useRef<HTMLCanvasElement>(null);

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    const ctx = canvas.getContext('2d');
    if (!ctx) return;

    // Draw saturation gradient (left to right)
    const saturationGradient = ctx.createLinearGradient(0, 0, canvas.width, 0);
    saturationGradient.addColorStop(0, `hsl(${hue}, 0%, 50%)`);
    saturationGradient.addColorStop(1, `hsl(${hue}, 100%, 50%)`);
    ctx.fillStyle = saturationGradient;
    ctx.fillRect(0, 0, canvas.width, canvas.height);

    // Draw lightness gradient (top to bottom)
    const lightnessGradient = ctx.createLinearGradient(0, 0, 0, canvas.height);
    lightnessGradient.addColorStop(0, 'rgba(255, 255, 255, 1)');
    lightnessGradient.addColorStop(0.5, 'rgba(255, 255, 255, 0)');
    lightnessGradient.addColorStop(0.5, 'rgba(0, 0, 0, 0)');
    lightnessGradient.addColorStop(1, 'rgba(0, 0, 0, 1)');
    ctx.fillStyle = lightnessGradient;
    ctx.fillRect(0, 0, canvas.width, canvas.height);
  }, [hue]);

  const { setSaturation, setLightness } = useColorPicker();

  const handleMove = useCallback(
    (e: React.MouseEvent<HTMLCanvasElement>) => {
      const canvas = canvasRef.current;
      if (!canvas) return;

      const rect = canvas.getBoundingClientRect();
      const x = Math.max(0, Math.min(e.clientX - rect.left, rect.width));
      const y = Math.max(0, Math.min(e.clientY - rect.top, rect.height));

      const saturationValue = (x / rect.width) * 100;
      const lightnessValue = 100 - (y / rect.height) * 100;

      setPos({ x: saturationValue, y: 100 - lightnessValue });
      setSaturation(saturationValue);
      setLightness(lightnessValue);
    },
    [setSaturation, setLightness]
  );

  const handleMouseDown = (e: React.MouseEvent<HTMLCanvasElement>) => {
    handleMove(e);
    const handleMouseMove = (e: MouseEvent) => {
      handleMove(e as any);
    };
    const handleMouseUp = () => {
      document.removeEventListener('mousemove', handleMouseMove);
      document.removeEventListener('mouseup', handleMouseUp);
    };
    document.addEventListener('mousemove', handleMouseMove);
    document.addEventListener('mouseup', handleMouseUp);
  };

  return (
    <div className="relative h-48 w-full cursor-crosshair rounded-md overflow-hidden">
      <canvas
        ref={canvasRef}
        width={300}
        height={192}
        className="h-full w-full"
        onMouseDown={handleMouseDown}
      />
      <div
        className="absolute h-4 w-4 -translate-x-1/2 -translate-y-1/2 rounded-full border-2 border-white shadow-md pointer-events-none"
        style={{
          left: `${pos.x}%`,
          top: `${pos.y}%`,
        }}
      />
    </div>
  );
};

export const ColorPickerHue = () => {
  const { hue, setHue } = useColorPicker();

  return (
    <Slider
      value={hue}
      onChange={(value) => setHue(value as number)}
      minValue={0}
      maxValue={360}
      step={1}
      aria-label="Hue"
      className="w-full"
      classNames={{
        track: "h-3 rounded-md",
        filler: "h-3 rounded-md",
        thumb: "h-5 w-5 bg-white border-2 border-gray-300"
      }}
      renderThumb={(props) => (
        <div
          {...props}
          className="h-5 w-5 rounded-full bg-white border-2 border-gray-300 shadow-md"
        />
      )}
      style={{
        background: 'linear-gradient(to right, #ff0000 0%, #ffff00 17%, #00ff00 33%, #00ffff 50%, #0000ff 67%, #ff00ff 83%, #ff0000 100%)',
      } as any}
    />
  );
};

export const ColorPickerAlpha = () => {
  const { alpha, setAlpha } = useColorPicker();

  return (
    <Slider
      value={alpha * 100}
      onChange={(value) => setAlpha((value as number) / 100)}
      minValue={0}
      maxValue={100}
      step={1}
      aria-label="Alpha"
      className="w-full"
      classNames={{
        track: "h-3 rounded-md",
        filler: "h-3 rounded-md",
        thumb: "h-5 w-5 bg-white border-2 border-gray-300"
      }}
    />
  );
};

export const ColorPickerOutput = () => {
  const { hue, saturation, lightness, mode } = useColorPicker();
  const color = useMemo(
    () => Color.hsl(hue, saturation, lightness),
    [hue, saturation, lightness]
  );

  const outputValue = useMemo(() => {
    switch (mode) {
      case 'hex':
        return color.hex(); // No alpha
      case 'rgb':
        return color.rgb().string();
      case 'hsl':
      default:
        return color.hsl().string();
    }
  }, [color, mode]);

  return (
    <Input
      value={outputValue}
      readOnly
      size="sm"
      classNames={{
        input: "text-xs font-mono"
      }}
    />
  );
};

export const ColorPickerFormat = () => {
  const { mode, setMode } = useColorPicker();

  return (
    <Select
      size="sm"
      selectedKeys={[mode]}
      onSelectionChange={(keys) => {
        const selected = Array.from(keys)[0];
        if (selected) setMode(selected.toString());
      }}
      className="w-24"
    >
      <SelectItem key="hsl">HSL</SelectItem>
      <SelectItem key="rgb">RGB</SelectItem>
      <SelectItem key="hex">HEX</SelectItem>
    </Select>
  );
};

export const ColorPickerEyeDropper = () => {
  const { setHue, setSaturation, setLightness, setAlpha } = useColorPicker();

  const handleEyeDropper = async () => {
    if (!('EyeDropper' in window)) {
      alert('EyeDropper API is not supported in your browser');
      return;
    }

    try {
      const eyeDropper = new (window as any).EyeDropper();
      const result = await eyeDropper.open();
      const color = Color(result.sRGBHex);
      setHue(color.hue());
      setSaturation(color.saturationl());
      setLightness(color.lightness());
      setAlpha(color.alpha());
    } catch (e) {
      // User cancelled
    }
  };

  return (
    <Button
      size="sm"
      variant="bordered"
      isIconOnly
      onPress={handleEyeDropper}
    >
      <Pipette className="h-4 w-4" />
    </Button>
  );
};
