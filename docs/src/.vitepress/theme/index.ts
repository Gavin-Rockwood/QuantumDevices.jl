import DefaultTheme from 'vitepress/theme'
import type { Theme } from 'vitepress'
import VersionPicker from '@/VersionPicker.vue'
import './style.css'
import './docstrings.css'

export default {
  extends: DefaultTheme,
  enhanceApp({ app }) { app.component('VersionPicker', VersionPicker) },
} satisfies Theme
