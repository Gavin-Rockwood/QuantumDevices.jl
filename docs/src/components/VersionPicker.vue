<script setup lang="ts">
import { computed, onMounted, ref } from 'vue'
import { useData, useRoute } from 'vitepress'
import VPNavScreenMenuGroup from 'vitepress/dist/client/theme-default/components/VPNavScreenMenuGroup.vue'
import VPNavBarMenuGroup from 'vitepress/dist/client/theme-default/components/VPNavBarMenuGroup.vue'

declare global {
  interface Window {
    DOC_VERSIONS?: string[]
    DOCUMENTER_CURRENT_VERSION?: string
  }
}
defineProps<{ screenMenu?: boolean }>()
const { site } = useData()
const route = useRoute()
const current = ref('dev')
const versions = ref<string[]>([])
const versionRoot = ref('')
const items = computed(() => {
  const path = route.path.startsWith(site.value.base)
    ? route.path.slice(site.value.base.length) : ''
  const suffix = typeof window === 'undefined' ? '' : window.location.search + window.location.hash
  return versions.value.length
    ? versions.value.map(version => ({
      text: version, link: `${versionRoot.value}/${version}/${path}${suffix}`,
      target: '_self', noIcon: true,
    }))
    : [{ text: 'Local preview', link: site.value.base, target: '_self', noIcon: true }]
})

async function switchVersion(event: MouseEvent) {
  const link = (event.target as Element).closest<HTMLAnchorElement>('a')
  if (!link || !versions.value.length || event.button !== 0 ||
      event.ctrlKey || event.metaKey || event.shiftKey || event.altKey) return
  event.preventDefault()
  const destination = new URL(link.href)
  const version = items.value.find(item => item.link === link.href)?.text
  if (!version) return
  // Probe the actual HTML file: older releases need no theme changes.
  const probe = new URL(destination)
  probe.search = ''
  probe.hash = ''
  if (probe.pathname.endsWith('/')) probe.pathname += 'index.html'
  else if (!probe.pathname.endsWith('.html')) probe.pathname += '.html'
  try {
    const response = await fetch(probe, { method: 'HEAD', signal: AbortSignal.timeout(5000) })
    if (response.status === 404 || response.status === 410) {
      const fallback = new URL(`${site.value.base}version-unavailable`, window.location.origin)
      fallback.searchParams.set('version', version)
      fallback.searchParams.set('page', destination.pathname)
      window.location.assign(fallback.href)
      return
    }
  } catch {
    // A network error is not evidence that the page is absent; try normal navigation.
  }
  window.location.assign(destination.href)
}
onMounted(() => {
  versionRoot.value = new URL('..', new URL(site.value.base, window.location.origin)).href.replace(/\/$/, '')
  current.value = window.DOCUMENTER_CURRENT_VERSION ?? 'dev'
  versions.value = window.DOC_VERSIONS ?? []
})
</script>

<template>
  <div @click.capture="switchVersion">
    <VPNavScreenMenuGroup v-if="screenMenu" :text="current" :items="items" />
    <VPNavBarMenuGroup v-else :item="{ text: current, items }" />
  </div>
</template>
