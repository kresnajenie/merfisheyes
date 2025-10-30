# How to Upload Chunked Datasets

## Step-by-Step Guide

### 1️⃣ Process Your H5AD File with Python

```bash
# In your terminal
python scripts/process_h5ad.py my_dataset.h5ad output_folder/
```

This creates a folder structure like:
```
output_folder/
├── manifest.json
├── coords/
│   ├── spatial.bin.gz
│   └── umap.bin.gz
├── expr/
│   ├── index.json
│   ├── chunk_00000.bin.gz
│   └── ...
├── obs/
│   └── ...
└── palettes/
    └── ...
```

### 2️⃣ Go to MERFISH Eyes Homepage

Open your browser and navigate to: `http://localhost:3000` (or your deployment URL)

### 3️⃣ Make Sure You're in "Single Cell" Mode

At the top of the page, you'll see a toggle switch:

```
[single cell] ⚪──────○ [single molecule]
                ^
              (should be here - blue/left side)
```

### 4️⃣ Click on the "Chunked Folder" Upload Box

You'll see 4 upload boxes:

```
┌───────────────┐  ┌───────────────┐  ┌───────────────┐  ┌───────────────┐
│  H5AD File    │  │Chunked Folder │  │ Xenium Folder │  │Merscope Folder│
│               │  │               │  │               │  │               │
│ Single .h5ad  │  │Pre-processed  │  │ Select Xenium │  │Select Merscope│
│     file      │  │chunked folder │  │output folder  │  │output folder  │
│               │  │               │  │               │  │               │
│               │  │ 👈 CLICK HERE │  │               │  │               │
└───────────────┘  └───────────────┘  └───────────────┘  └───────────────┘
```

**Important:** Use "Chunked Folder" for your Python-processed datasets, NOT "H5AD File"!

### 5️⃣ Select Your Chunked Folder

When the file picker opens:

**For macOS/Linux:**
- Navigate to where you saved `output_folder/`
- **Select the entire folder** (not the files inside)
- Click "Upload" or "Open"

**For Windows:**
- Navigate to where you saved `output_folder/`
- **Select the folder**
- Click "Select Folder" or "Upload"

**Important:** You need to select the **entire folder**, not individual files!

### 6️⃣ Automatic Detection

The browser will automatically:
1. ✅ Detect it's a chunked dataset (by checking for `manifest.json`, `expr/index.json`, etc.)
2. ✅ Load it using `ChunkedDataAdapter` in local mode
3. ✅ Show progress: "Detected chunked dataset folder..."
4. ✅ Navigate to `/viewer` when complete

### 7️⃣ Start Visualizing!

Once loaded, you can:
- Select genes from the gene panel (loads expression on-demand from chunks)
- Filter by cell types
- Explore embeddings (UMAP, etc.)
- Adjust visualization settings

---

## 🎯 What Files Are Detected as Chunked?

The browser checks for these required files:
- ✅ `manifest.json`
- ✅ `expr/index.json`
- ✅ `expr/chunk_*.bin.gz` (at least one chunk file)
- ✅ `coords/spatial.bin.gz`

If all are present → **Chunked dataset detected** → Instant loading!

If not → Falls back to normal H5AD file processing (if you selected a `.h5ad` file)

---

## ⚡ Benefits of Chunked Upload

| Feature | Regular H5AD Upload | Chunked Folder Upload |
|---------|-------------------|----------------------|
| **Browser Processing** | ✅ Full parsing in browser | ❌ None (pre-processed) |
| **Load Time** | 10-60 seconds | 1-3 seconds |
| **Memory Usage** | Full matrix in RAM | Metadata only |
| **Gene Expression** | Pre-loaded | On-demand from chunks |
| **Dataset Size Limit** | ~500K cells (browser limit) | ♾️ Unlimited |

---

## 🔍 Troubleshooting

### "No files detected" or "Invalid dataset"
- Make sure you selected the **entire folder**, not files inside
- Check that all required files are present (manifest.json, expr/index.json, etc.)
- Verify the Python script completed successfully

### "Failed to load dataset"
- Check browser console for errors (F12 → Console tab)
- Make sure the folder structure matches the expected format
- Try re-running the Python script

### Folder selection not working
- **Chrome/Edge**: Full folder upload support ✅
- **Firefox**: Full folder upload support ✅
- **Safari**: May have limited folder support - try Chrome/Edge

### File picker only shows files, not folders
- Look for a "Select Folder" or "Upload Folder" option in the file picker
- Some browsers require you to open the folder first, then click "Upload" on the folder itself

---

## 📝 Example Workflow

```bash
# 1. Process your dataset
cd /path/to/merfisheyes-heroui
python scripts/process_h5ad.py ~/data/brain_large.h5ad ~/data/brain_processed/

# 2. Start the app (if not running)
npm run dev

# 3. Open browser
# http://localhost:3000

# 4. Upload
# Click "H5AD" box → Select ~/data/brain_processed/ folder

# 5. Visualize
# Browser navigates to /viewer automatically
# Select genes, explore data!
```

---

## 🎉 That's It!

You can now process any size H5AD dataset with Python and load it instantly in the browser without memory limitations!
