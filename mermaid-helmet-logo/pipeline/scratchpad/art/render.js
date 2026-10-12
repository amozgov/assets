// Render SVG files to PNG next to them (same basename).
const { chromium } = require('playwright');
const fs = require('fs');
const path = require('path');
(async () => {
  const files = process.argv.slice(2);
  const browser = await chromium.launch();
  const page = await browser.newPage({ deviceScaleFactor: 1 });
  for (const f of files) {
    const svg = fs.readFileSync(f, 'utf8');
    const w = +(svg.match(/width="([\d.]+)"/) || [0, 800])[1];
    const h = +(svg.match(/height="([\d.]+)"/) || [0, 800])[1];
    await page.setViewportSize({ width: Math.ceil(w), height: Math.ceil(h) });
    await page.goto('file://' + path.resolve(f));
    await page.screenshot({ path: f.replace(/\.svg$/, '.png'), omitBackground: false });
  }
  await browser.close();
})();
