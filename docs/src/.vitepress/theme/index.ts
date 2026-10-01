import DefaultTheme from 'vitepress/theme'
import type { Theme } from 'vitepress'
import VersionPicker from '@/VersionPicker.vue'
import VersionUnavailable from '@/VersionUnavailable.vue'
import './style.css'
import './docstrings.css'

export default {
  extends: DefaultTheme,
  enhanceApp({ app }) {
    app.component('VersionPicker', VersionPicker)
    app.component('VersionUnavailable', VersionUnavailable)
  },
} satisfies Theme
