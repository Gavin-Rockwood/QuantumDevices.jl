<script setup lang="ts">
import { onMounted, ref } from 'vue'
import { useData } from 'vitepress'

const { site } = useData()
const version = ref('the selected version')
const home = ref(site.value.base)
onMounted(() => {
  const selected = new URLSearchParams(window.location.search).get('version')
  if (selected && /^(?:v?\d+\.\d+(?:\.\d+)?(?:[-+][\w.-]+)?|dev|stable|previews\/PR\d+)$/.test(selected)) {
    version.value = selected
    home.value = new URL(`${encodeURIComponent(selected)}/`,
      new URL('..', new URL(site.value.base, window.location.origin))).href
  }
})
</script>

<template>
  <p>This page doesn't exist in {{ version }}.</p>
  <p><a :href="home" target="_self">Take me home</a></p>
</template>
