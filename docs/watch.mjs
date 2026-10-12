import { watch } from 'node:fs'
import { spawn } from 'node:child_process'
import path from 'node:path'
import { fileURLToPath } from 'node:url'

const docs = path.dirname(fileURLToPath(import.meta.url))
const root = path.dirname(docs)
// Literate writes these files during generation; their Julia sources are watched.
const generated = new Set([
  'tutorials/transmon_resonator_control.md',
  'tutorials/tunable_coupler_control.md',
])
let generator, server, timer
let pending = false
let stopping = false

function startServer() {
  const child = spawn(process.execPath, [
    path.join(docs, 'node_modules/vitepress/bin/vitepress.js'),
    'dev', 'build/.documenter', '--host', '127.0.0.1',
  ], { cwd: docs, stdio: 'inherit' })
  server = child
  child.on('error', error => { console.error(error.message); stop(1) })
  child.on('exit', code => { if (!stopping && server === child) stop(code ?? 1) })
}

function regenerate() {
  if (stopping || generator) return
  pending = false
  // Documenter rewrites VitePress configuration in several steps. Restart only
  // after generation finishes so VitePress never reads a partial configuration.
  const previousServer = server
  server = undefined
  previousServer?.kill('SIGTERM')
  console.log('\nRegenerating docs (examples and doctests skipped)…')
  generator = spawn('julia', ['--startup-file=no', '--project=docs', 'docs/make.jl'], {
    cwd: root,
    stdio: 'inherit',
    env: { ...process.env, DOCS_DRAFT: 'true', DOCS_MARKDOWN_ONLY: 'true', DOCS_DEPLOY: 'false' },
  })
  generator.on('error', error => { console.error(error.message); stop(1) })
  generator.on('exit', code => {
    generator = undefined
    if (stopping) return
    if (code === 0 && !server) startServer()
    if (code !== 0) console.error('Generation failed. Fix the source and save to retry.')
    if (pending) regenerate()
  })
}

function schedule() {
  pending = true
  clearTimeout(timer)
  timer = setTimeout(regenerate, 300)
}

const watchers = [
  watch(path.join(docs, 'src'), { recursive: true }, (_event, filename) => {
    // Directory events also fire when Literate writes its generated pages.
    if (filename && path.extname(filename) && !generated.has(filename.replaceAll(path.sep, '/'))) schedule()
  }),
  watch(path.join(docs, 'make.jl'), schedule),
]
for (const watcher of watchers) {
  watcher.on('error', error => { console.error(`File watching failed: ${error.message}`); stop(1) })
}

function stop(code = 0) {
  if (stopping) return
  stopping = true
  clearTimeout(timer)
  for (const watcher of watchers) watcher.close()
  generator?.kill('SIGTERM')
  server?.kill('SIGTERM')
  process.exitCode = code
}
process.on('SIGINT', () => stop())
process.on('SIGTERM', () => stop())
regenerate()
