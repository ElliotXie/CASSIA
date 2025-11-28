import fs from "fs"
import path from "path"
import matter from "gray-matter"

const vignetteDirectory = path.join(process.cwd(), "content/vignette")

// Create the directory if it doesn't exist
try {
  if (!fs.existsSync(vignetteDirectory)) {
    fs.mkdirSync(vignetteDirectory, { recursive: true });
  }
} catch (error) {
  console.error("Error creating vignette directory:", error);
}

export function getVignetteBySlug(slug: string) {
  const realSlug = slug.replace(/\.md$/, "")
  const fullPath = path.join(vignetteDirectory, `${realSlug}.md`)

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
    }
  } catch (error) {
    console.error(`Error reading file ${fullPath}:`, error)
    return {
      slug: realSlug,
      frontmatter: { title: "Not Found" },
      content: "# Not Found\n\nThe requested vignette could not be found.",
    }
  }
}

export function getAllVignettes() {
  try {
    if (!fs.existsSync(vignetteDirectory)) {
      return [];
    }
    
    const slugs = fs
      .readdirSync(vignetteDirectory)
      .filter((file) => file.endsWith(".md"))
      .map((file) => file.replace(/\.md$/, ""))

    const vignettes = slugs.map((slug) => getVignetteBySlug(slug))

    return vignettes
  } catch (error) {
    console.error("Error reading vignette directory:", error)
    return []
  }
} 