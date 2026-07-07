import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

// base './' + hash routing => works on GitHub Pages project sites without
// any path rewriting (visualization_plan.md §7 decision 4).
export default defineConfig({
  base: './',
  plugins: [react()],
})
