import fs from "fs"
import path from "path"
import matter from "gray-matter"

const contentDirectory = path.join(process.cwd(), "content")

export function getVignetteBySlug(slug: string, locale: string = 'en') {
  const realSlug = slug.replace(/\.md$/, "")

  // Try locale-specific path first
  let fullPath = path.join(contentDirectory, locale, "vignette", `${realSlug}.md`)
  let isFallback = false

  // Fallback to English if locale-specific file doesn't exist
  if (!fs.existsSync(fullPath) && locale !== 'en') {
    fullPath = path.join(contentDirectory, 'en', "vignette", `${realSlug}.md`)
    isFallback = true
  }

  try {
    const fileContents = fs.readFileSync(fullPath, "utf8")
    const { data, content } = matter(fileContents)

    // Fix the triple backtick escaping in markdown
    const processedContent = content
      // Replace escaped backticks with actual backticks
      .replace(/\\`\\`\\`/g, "```")
      // Make sure only content inside triple backticks is treated as code blocks
      .replace(/```(\s*)(r|R)?\s*\n([\s\S]*?)```/g, (match, space, lang, code) => {
        return `\`\`\`${lang || "r"}\n${code}\`\`\``
      })

    return {
      slug: realSlug,
      frontmatter: data,
      content: processedContent,
      isFallback,
    }
  } catch (error) {
    console.error(`Error reading file ${fullPath}:`, error)
    return {
      slug: realSlug,
      frontmatter: { title: "Not Found" },
      content: "# Not Found\n\nThe requested vignette could not be found.",
      isFallback: false,
    }
  }
}

export function getAllVignettes(locale: string = 'en') {
  // Get vignettes for the locale, with English fallback
  const enVignetteDir = path.join(contentDirectory, 'en', 'vignette')

  try {
    if (!fs.existsSync(enVignetteDir)) {
      return [];
    }

    // Use English slugs as the canonical list
    const slugs = fs
      .readdirSync(enVignetteDir)
      .filter((file) => file.endsWith(".md"))
      .map((file) => file.replace(/\.md$/, ""))

    const vignettes = slugs.map((slug) => getVignetteBySlug(slug, locale))

    return vignettes
  } catch (error) {
    console.error("Error reading vignette directory:", error)
    return []
  }
} 