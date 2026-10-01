import { defineConfig } from 'vitepress'
import path from 'node:path'

// Documenter replaces these placeholders for branch, release, and preview builds.
const location = { base: 'REPLACE_ME_DOCUMENTER_VITEPRESS' }
const repositoryBase = '/' + location.base.split('/').filter(Boolean)[0] + '/'
const versionScripts: any[] = process.env.DOCS_DEPLOY === 'true' ? [
  ['script', { src: `${repositoryBase}versions.js` }],
  ['script', { src: `${location.base}siteinfo.js` }],
] : []

export default defineConfig({
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  title: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  description: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  outDir: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  // Documenter parses Markdown before VitePress; supply homepage metadata here.
  transformPageData(page) {
    if (page.relativePath === 'version-unavailable.md') {
      page.frontmatter = { ...page.frontmatter, sidebar: false, aside: false, prev: false, next: false }
      return
    }
    if (page.relativePath !== 'index.md') return
    page.frontmatter = {
      ...page.frontmatter,
      layout: 'home',
      hero: {
        name: 'QuantumDevices.jl',
        text: 'From device Hamiltonians to calibrated gates',
        tagline: 'Build symbolic quantum-device models, retain the physics you need, and simulate controlled evolution in Julia.',
        actions: [
          { theme: 'brand', text: 'Start building', link: '/getting_started/quickstart' },
          { theme: 'alt', text: 'API reference', link: '/resources/api' },
        ],
      },
      features: [
        { title: 'Compose device models', details: 'Combine qubits and transmons using symbolic parameters and ordered operator products.' },
        { title: 'Shape and calibrate controls', details: 'Use built-in pulses, independent flattop ramps, and caller-selected SciML optimizers.' },
        { title: 'Keep results reproducible', details: 'Track dressed states and save reconstructible models with their named gates.' },
      ],
    }
  },
  lastUpdated: true,
  cleanUrls: true,
  head: [
    ['link', { rel: 'icon', type: 'image/svg+xml', href: `${location.base}favicon.svg` }],
    ...versionScripts,
  ],
  vite: { resolve: { alias: { '@': path.resolve(__dirname, '../components') } } },
  markdown: {
    math: true,
    theme: { light: 'github-light', dark: 'github-dark' },
  },
  themeConfig: {
    outline: [2, 3],
    search: { provider: 'local', options: { detailedView: true } },
    nav: [
      { text: 'Start', link: '/getting_started/quickstart' },
      { text: 'Guide', link: '/user_guide/models' },
      {
        text: 'Tutorials',
        activeMatch: '^/tutorials/',
        items: [
          { text: 'Transmon resonator control', link: '/tutorials/transmon_resonator_control' },
          { text: 'Tunable coupler control', link: '/tutorials/tunable_coupler_control' },
        ],
      },
      { text: 'API', link: '/resources/api' },
      { text: 'Development', link: '/development/extensions' },
      { component: 'VersionPicker' },
    ],
    sidebar: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    editLink: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    socialLinks: [{ icon: 'github', link: 'REPLACE_ME_DOCUMENTER_VITEPRESS' }],
    footer: { message: 'Built with Julia, Documenter, and VitePress.' },
  },
})
