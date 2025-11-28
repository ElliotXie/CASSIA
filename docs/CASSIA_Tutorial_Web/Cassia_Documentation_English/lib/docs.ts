import fs from "fs"
import path from "path"
import matter from "gray-matter"

const contentDirectory = path.join(process.cwd(), "content")

export function getDocBySlug(slug: string, locale: string = 'en') {
  const realSlug = slug.replace(/\.md$/, "")

  // Try locale-specific path first
  let fullPath = path.join(contentDirectory, locale, "docs", `${realSlug}.md`)
  let isFallback = false

  // Fallback to English if locale-specific file doesn't exist
  if (!fs.existsSync(fullPath) && locale !== 'en') {
    fullPath = path.join(contentDirectory, 'en', "docs", `${realSlug}.md`)
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
      content: "# Not Found\n\nThe requested document could not be found.",
      isFallback: false,
    }
  }
}

export function getAllDocs(locale: string = 'en') {
  // Get docs for the locale, with English fallback
  const enDocsDir = path.join(contentDirectory, 'en', 'docs')

  try {
    // Use English slugs as the canonical list
    const slugs = fs
      .readdirSync(enDocsDir)
      .filter((file) => file.endsWith(".md"))
      .map((file) => file.replace(/\.md$/, ""))

    const docs = slugs.map((slug) => getDocBySlug(slug, locale))

    return docs
  } catch (error) {
    console.error("Error reading docs directory:", error)
    return []
  }
}
