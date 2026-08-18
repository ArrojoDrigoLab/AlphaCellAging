// ================================================
// QuPath Groovy script:
// Export selected measurements (including distances) as CSV per image
// ================================================

import qupath.lib.gui.scripting.QPEx

// ——— USER PARAMETERS ———
def desiredColumns = [
    'Image','Object ID','Object type','Name','Classification',
    'Centroid X µm','Centroid Y µm',
    'DAPI: Cell: Mean','DAPI: Cell: Median',
    'GCG: Cell: Mean','GCG: Cell: Median',
    'CD45: Cell: Mean','CD45: Cell: Median',
    'CD20: Cell: Mean','CD20: Cell: Median',
    'CD19: Cell: Mean','CD19: Cell: Median',
    'IRX2: Cell: Mean','IRX2: Nucleus: Mean',
    'CD8a: Cell: Mean','CD8a: Cell: Median',
    'CD8: Cell: Mean','CD8: Cell: Median',
    'LAM: Cell: Mean','LAM: Nucleus: Mean',
    'PDX1: Cell: Mean','PDX1: Nucleus: Median',
    'NEUROD1: Cell: Mean','NEUROD1: Nucleus: Mean',
    'CD3e: Cell: Mean','CD3e: Cell: Median',
    'CD3: Cell: Mean','CD3: Cell: Median',
    'SST: Cell: Mean','SST: Cell: Median',
    'ECAD: Cell: Mean','ECAD: Cell: Median',
    'HLA-DR: Cell: Mean','HLA-DR: Cell: Median',
    'LAMA2: Cell: Mean','LAMA2: Cell: Median',
    'CD4: Cell: Mean','CD4: Cell: Median',
    'ACTA2: Cell: Mean','ACTA2: Cell: Median',
    'COL1A1: Cell: Mean','COL1A1: Cell: Median',
    'VIM: Cell: Mean','VIM: Cell: Median',
    'proINS: Cell: Mean','proINS: Cell: Median',
    'PAX6: Cell: Mean','PAX6: Nucleus: Mean',
    'PPY: Cell: Mean','PPY: Cell: Median',
    'CD11c: Cell: Mean','CD11c: Cell: Median',
    'NKX6-1: Cell: Mean','NKX6-1: Nucleus: Mean',
    'KRT: Cell: Mean','KRT: Cell: Median',
    'ARX: Cell: Mean','ARX: Nucleus: Mean',
    'EPCAM: Cell: Mean','EPCAM: Cell: Median',
    'SELP: Cell: Mean','SELP: Cell: Median',
    'CD31: Cell: Mean','CD31: Cell: Median',
    'HSPG2: Cell: Mean','HSPG2: Cell: Median',
    'CD117: Cell: Mean','CD117: Cell: Median',
    'NKX2-2: Cell: Mean','NKX2-2: Nucleus: Mean',
    'GP2: Cell: Mean','GP2: Cell: Median',
    'CDX2: Cell: Mean','CDX2: Nucleus: Mean',
    'GHRL: Cell: Mean','GHRL: Cell: Median',
    'MMR: Cell: Mean','MMR: Cell: Median',
    'ATP1A1: Cell: Mean','ATP1A1: Cell: Median',
    'FN1: Cell: Mean','FN1: Cell: Median',
    'PD-L1: Cell: Mean','PD-L1: Cell: Median',
    'MPO: Cell: Mean','MPO: Cell: Median',
    'IBA1: Cell: Mean','IBA1: Cell: Median',
    'MCAM: Cell: Mean','MCAM: Cell: Median',
    'CHGA: Cell: Mean','CHGA: Cell: Median',
    'SYP: Cell: Mean','SYP: Cell: Median',
    'CD163: Cell: Mean','CD163: Cell: Median',
    'CD68: Cell: Mean','CD68: Cell: Median',
    'COL6: Cell: Mean','COL6: Cell: Median',
    'NPY: Cell: Mean','NPY: Cell: Median',
    'COL4A1: Cell: Mean','COL4A1: Cell: Median',
    'CD66: Cell: Mean','CD66: Cell: Median',
    'CD66b: Cell: Mean','CD66b: Cell: Median',
    'CD44: Cell: Mean','CD44: Cell: Median',
    'PNLIP: Cell: Mean','PNLIP: Cell: Median',
    'CD141: Cell: Mean','CD141: Cell: Median',
    'CD90: Cell: Mean','CD90: Cell: Median',
    'CD39L3: Cell: Mean','CD39L3: Cell: Median',
    'SOX9: Cell: Mean','SOX9: Nucleus: Mean',
    'CPEP: Cell: Mean','CPEP: Cell: Median',
    'Ki67: Cell: Mean','Ki67: Nucleus: Mean',
    'TUBB3: Cell: Mean','TUBB3: Cell: Median',
    'CD40: Cell: Mean','CD40: Cell: Median',
    'LYVE1: Cell: Mean','LYVE1: Cell: Median',
    'DPP4: Cell: Mean','DPP4: Cell: Median'
]

// — Specify an absolute export folder here —
def exportDir = new File('D:/Mike/CODEX/T2DM Donors')
if (!exportDir.exists()) exportDir.mkdirs()

// Get current image and name
def imageData  = QPEx.getCurrentImageData()
def imageName  = imageData.getServer().getMetadata().getName()

// Build output CSV file path
def outFile = new File(exportDir, imageName + '_measurements.csv')

// 1) Write header row
outFile.withPrintWriter { pw ->
    pw.println(desiredColumns.join(','))
}

// 2) Collect all detections (swap to getAnnotationObjects() if needed)
def objects = getDetectionObjects()
if (objects.isEmpty()) {
    println "No detections for image ${imageName}; skipping."
    return
}

// 3) Write one CSV row per object
objects.eachWithIndex { obj, idx ->
    def ml = obj.getMeasurementList()
    def roi = obj.getROI()
    def row = desiredColumns.collect { col ->
        switch (col) {
            case 'Image':
                return imageName
            case 'Object ID':
                return (idx + 1).toString()
            case 'Object type':
                return obj.getClass().getSimpleName()
            case 'Name':
                return obj.getName() ?: ''
            case 'Classification':
                return obj.getPathClass()?.toString() ?: ''
            case ['Centroid X µm', 'Centroid X um']:
                return roi != null ? roi.getCentroidX() : ''
            case ['Centroid Y µm', 'Centroid Y um']:
                return roi != null ? roi.getCentroidY() : ''
            default:
                def raw = ml.get(col)
                return raw != null ? raw.toString() : ''
        }
    }.join(',')
    outFile.append(row + '\n')
}

// 4) Final confirmation
println "Wrote ${objects.size()} rows × ${desiredColumns.size()} columns to ${outFile}"
