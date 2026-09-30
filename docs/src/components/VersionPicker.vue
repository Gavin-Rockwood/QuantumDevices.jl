<script setup lang="ts">
import { computed, onMounted, ref } from 'vue'
import { useData } from 'vitepress'
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
const current = ref('dev')
const versions = ref<string[]>([])
const versionRoot = ref('')
const items = computed(() => versions.value.length
  ? versions.value.map(version => ({ text: version, link: `${versionRoot.value}/${version}/` }))
  : [{ text: 'Local preview', link: '/' }])
onMounted(() => {
  versionRoot.value = `${window.location.origin}/${site.value.base.split('/').filter(Boolean)[0]}`
  current.value = window.DOCUMENTER_CURRENT_VERSION ?? 'dev'
  versions.value = window.DOC_VERSIONS ?? []
})
</script>

<template>
  <VPNavScreenMenuGroup v-if="screenMenu" :text="current" :items="items" />
  <VPNavBarMenuGroup v-else :item="{ text: current, items }" />
</template>
