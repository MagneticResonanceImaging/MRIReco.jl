import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import { mathjaxPlugin } from './mathjax-plugin'
import { juliaReplTransformer } from './julia-repl-transformer'
import footnote from "markdown-it-footnote";
import path from 'path'

const mathjax = mathjaxPlugin()

function getBaseRepository(base: string): string {
  if (!base || base === '/') return '/';
  const parts = base.split('/').filter(Boolean);
  return parts.length > 0 ? `/${parts[0]}/` : '/';
}

const baseTemp = {
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',// TODO: replace this in makedocs!
}

const navTemp = {
  nav: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
}

// The sections of `pages` in make.jl. The site title links to the home page, so the top bar leaves it out.
const sections = navTemp.nav.filter((item: any) => item.link !== '/index')

const nav = [
  ...sections,
  {
    component: 'VersionPicker'
  }
]

// The links of the pages of a section
const links = (item: any): string[] => (item.link ? [item.link] : (item.items ?? []).flatMap(links))

// Each page shows the sidebar of its section, chosen by the folder of the page or, for a page at the top level of
// docs/src, by the page itself
const sidebar: Record<string, any[]> = {}
for (const section of sections) {
  for (const link of links(section)) {
    const folder = link.slice(0, link.lastIndexOf('/') + 1)
    const key = folder === '/' ? `${link}.md` : folder
    if (sidebar[key] && sidebar[key][0] !== section) throw new Error(`${key} holds pages of two sections`)
    sidebar[key] = [section]
  }
}

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',// TODO: replace this in makedocs!
  title: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  description: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  lastUpdated: true,
  cleanUrls: true,
  outDir: 'REPLACE_ME_DOCUMENTER_VITEPRESS', // This is required for MarkdownVitepress to work correctly...
  head: [
    ['link', { rel: 'icon', href: 'REPLACE_ME_DOCUMENTER_VITEPRESS_FAVICON' }],
    ['script', {src: `${getBaseRepository(baseTemp.base)}versions.js`}],
    // ['script', {src: '/versions.js'], for custom domains, I guess if deploy_url is available.
    ['script', {src: `${baseTemp.base}siteinfo.js`}],
    // REPLACE_ME_DOCUMENTER_VITEPRESS_NOINDEX
  ],

  markdown: {
    component: { inlineTags: ['CitationPreview'] },
    codeTransformers: [juliaReplTransformer()],
    config(md) {
      md.use(tabsMarkdownPlugin);
      md.use(footnote);
      mathjax.markdownConfig(md);
    },
    theme: {
      light: "github-light",
      dark: "github-dark"
    },
  },
  vite: {
    plugins: [
      mathjax.vitePlugin,
    ],
    define: {
      __DEPLOY_ABSPATH__: JSON.stringify('REPLACE_ME_DOCUMENTER_VITEPRESS_DEPLOY_ABSPATH'),
    },
    resolve: {
      alias: {
        '@': path.resolve(__dirname, '../components')
      }
    },
    optimizeDeps: {
      exclude: [
        '@nolebase/vitepress-plugin-enhanced-readabilities/client',
        'vitepress',
        '@nolebase/ui',
      ],
    },
    ssr: {
      noExternal: [
        // If there are other packages that need to be processed by Vite, you can add them here.
        '@nolebase/vitepress-plugin-enhanced-readabilities',
        '@nolebase/ui',
      ],
    },
  },
  themeConfig: {
    outline: 'deep',
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },
    nav,
    sidebar,
    sidebarDrawer: 'REPLACE_ME_DOCUMENTER_VITEPRESS_SIDEBAR_DRAWER',
    editLink: {
      // Literate.jl generates the pages in docs/src/generated from docs/src/literate, and make.jl copies the
      // contributing guide from CONTRIBUTING.md
      pattern: ({ filePath }: { filePath: string }): string => {
        const base = 'https://github.com/MagneticResonanceImaging/MRIReco.jl/edit/master'
        const generated = filePath.match(/^generated\/(.+)\.md$/)
        if (generated) return `${base}/docs/src/literate/${generated[1]}.jl`
        if (filePath === 'contributing.md') return `${base}/CONTRIBUTING.md`
        return `${base}/docs/src/${filePath}`
      }
    },
    socialLinks: [
      { icon: 'github', link: 'REPLACE_ME_DOCUMENTER_VITEPRESS' }
    ],
    footer: {
      message:
        'Made with <a href="https://documenter.juliadocs.org/stable/" target="_blank"><strong>Documenter.jl</strong></a>, ' +
        '<a href="https://vitepress.dev" target="_blank"><strong>VitePress</strong></a> and ' +
        '<a href="https://luxdl.github.io/DocumenterVitepress.jl/stable/" target="_blank"><strong>DocumenterVitepress.jl</strong></a><br>' +
        'Released under the MIT License.'
    }
  }
})
