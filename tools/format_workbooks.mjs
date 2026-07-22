import fs from "node:fs/promises";
import path from "node:path";
import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";

const root = "C:/Users/32574/Documents/Codex/2026-07-20/new-chat/ich-hfcg-reproducibility";
const targets = [
  path.join(root, "results/tables/TableS2_S4_raw_GEO_reanalysis.xlsx"),
  path.join(root, "results/tables/TableS5_S6_MR_colocalisation_complete.xlsx"),
];
const previewDir = path.join(root, "../work/spreadsheet_previews");
await fs.mkdir(previewDir, { recursive: true });

for (const target of targets) {
  const input = await FileBlob.load(target);
  const workbook = await SpreadsheetFile.importXlsx(input);
  const overview = await workbook.inspect({
    kind: "sheet",
    include: "id,name",
    maxChars: 8000,
  });
  console.log(`WORKBOOK ${path.basename(target)}`);
  console.log(overview.ndjson);
  const sheetRows = overview.ndjson.trim().split("\n").filter(Boolean).map(JSON.parse);
  for (const info of sheetRows) {
    const sheet = workbook.worksheets.getItemAt(info.index);
    const used = sheet.getUsedRange();
    sheet.showGridLines = false;
    sheet.freezePanes.freezeRows(1);
    used.format.font = { name: "Aptos", size: 10, color: "#1F2937" };
    used.format.columnWidth = 15;
    used.format.rowHeight = 18;
    const header = sheet.getRange(info.range.replace(/\d+$/, "1"));
    header.format = {
      fill: "#244A73",
      font: { name: "Aptos", size: 10, bold: true, color: "#FFFFFF" },
      verticalAlignment: "center",
      wrapText: true,
      borders: { bottom: { style: "medium", color: "#17324D" } },
    };
    header.format.rowHeight = 30;
    const maxPreviewRow = Math.min(15, Number(info.range.match(/(\d+)$/)[1]));
    const endColumn = info.range.match(/:([A-Z]+)/)[1];
    if (info.name.includes("Provenance") || info.name.includes("parameters")) {
      sheet.getRange(`B1:B${maxPreviewRow}`).format.columnWidth = 65;
      sheet.getRange(`B1:B${maxPreviewRow}`).format.wrapText = true;
    }
    const preview = await workbook.render({
      sheetName: info.name,
      range: `A1:${endColumn}${maxPreviewRow}`,
      scale: 1.2,
      format: "png",
    });
    await fs.writeFile(
      path.join(previewDir, `${path.basename(target, ".xlsx")}_${String(info.index + 1).padStart(2, "0")}.png`),
      new Uint8Array(await preview.arrayBuffer()),
    );
  }
  const errors = await workbook.inspect({
    kind: "match",
    searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
    options: { useRegex: true, maxResults: 100 },
    summary: "formula error scan",
  });
  console.log(errors.ndjson || "No formula errors");
  const output = await SpreadsheetFile.exportXlsx(workbook);
  const temp = `${target}.formatted.xlsx`;
  await output.save(temp);
  await fs.copyFile(temp, target);
  await fs.unlink(temp);
}
