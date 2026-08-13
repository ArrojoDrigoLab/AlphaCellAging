//————————————————————————————————————————————————————————
// Export only selected cell features per image — Run for project
//————————————————————————————————————————————————————————

import static qupath.lib.gui.scripting.QPEx.*       // for project & image context
import qupath.lib.projects.Project
import java.nio.charset.StandardCharsets

// 1) List exactly the measurements you want (as shown in Show Info…)
def keys = [
    'Centroid X px', 'Centroid Y px',
    'CD20: Cell: Mean',  'CD8: Cell: Mean',   'CD99: Cell: Mean',
    'NFkb: Cell: Mean',  'GnzB: Cell: Mean',  'Ki67: Cell: Mean',
    'PDX-1: Cell: Mean', 'CD56: Cell: Mean',  'Foxp3: Cell: Mean',
    'CD4: Cell: Mean',   'NKX6.1: Cell: Mean','Somatostatin: Cell: Mean',
    'CA2: Cell: Mean',   'PP: Cell: Mean',    'CD57: Cell: Mean',
    'CD31: Cell: Mean',  'CD14: Cell: Mean',  'C-peptide: Cell: Mean',
    'Nestin: Cell: Mean','pan-Keratin: Cell: Mean','CD44: Cell: Mean',
    'HLA-ABC: Cell: Mean','Glucagon: Cell: Mean','CD11b: Cell: Mean',
    'CD45: Cell: Mean','beta-Actin: Cell: Mean','CD68: Cell: Mean',
    'Collagen: Cell: Mean','CD3: Cell: Mean','pS6: Cell: Mean',
    'CD45RO: Cell: Mean','HLA-DR: Cell: Mean','GHRL: Cell: Mean'
]

// 2) Prepare output folder under the project directory
def projectDir = getProject().getProjectDirectory()
def outputDir  = new File(projectDir, 'CellMeasurementsFiltered')
outputDir.mkdirs()  // make it if it doesn't exist

// 3) This script runs **once per image** in the project
def imageData = getCurrentImageData()
def imageName = getImageName()
// Sanitize imageName so it makes a valid filename
def safeName  = imageName.replaceAll(/[\\/:*?"<>| ]+/, '_')

// Build the CSV file
def csvFile = new File(outputDir, safeName + '.csv')
csvFile.withPrintWriter(StandardCharsets.UTF_8.name()) { writer ->
    // Write header
    writer.println((['Image'] + keys).join(','))

    // Iterate each detection (cell)
    getDetectionObjects().each { cell ->
        def ml = cell.getMeasurementList()
        // Start row with image name
        def row = [imageName]
        // Append each requested measurement (or blank if missing)
        keys.each { key ->
            def v = ml.getMeasurement(key)
            row << (v != null ? v : '')
        }
        writer.println(row.join(','))
    }
}

print "✅ Exported ${imageName} → ${csvFile.absolutePath}"
