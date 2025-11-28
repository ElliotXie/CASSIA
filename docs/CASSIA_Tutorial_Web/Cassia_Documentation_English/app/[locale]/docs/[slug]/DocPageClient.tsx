"use client"

import ReactMarkdown, { Components } from "react-markdown"
import rehypeSlug from "rehype-slug"
import rehypeAutolinkHeadings from "rehype-autolink-headings"
import { useState } from "react"
import { notFound } from "next/navigation"
import { ClassAttributes, HTMLAttributes } from "react"

import { useTranslations } from "next-intl"

// Assuming DocType is the type for the doc object, define or import it
// If not defined elsewhere, create a basic type
type DocType = {
  slug: string
  frontmatter: { [key: string]: any; title: string }
  content: string
  isFallback?: boolean
}

// Helper function to create slug IDs for headings
function slugify(text: string) {
  return text
    .toString()
    .toLowerCase()
    .replace(/\s+/g, "-") // Replace spaces with -
    .replace(/[^\w-]+/g, "") // Remove all non-word chars
    .replace(/--+/g, "-") // Replace multiple - with single -
    .replace(/^-+/, "") // Trim - from start of text
    .replace(/-+$/, "") // Trim - from end of text
}

// Update the props interface to include 'doc'
interface DocPageClientProps {
  params: {
    slug: string
    locale: string
  }
  doc: DocType // Add the doc prop
}

// Define type for code component props provided by ReactMarkdown
interface CodeProps extends ClassAttributes<HTMLElement>, HTMLAttributes<HTMLElement> {
  node?: any; // You might want to type this more accurately if needed
  inline?: boolean;
  className?: string;
  children?: React.ReactNode;
}

// Function to slugify a string
// function slugify(string: string | number | boolean) {
//   return String(string)
//     .normalize('NFKD') // Remove diacritics
//     .toLowerCase()
//     .trim()
//     .replace(/\s+/g, '-') // Replace spaces with hyphens
//     .replace(/[^\w-]+/g, '') // Remove all non-word characters
//     .replace(/--+/g, '-') // Replace multiple hyphens with single hyphen
//     .replace(/^-+/, '') // Trim - from start of text
//     .replace(/-+$/, ''); // Trim - from end of text
// }

export default function DocPageClient({ params, doc }: DocPageClientProps) {
  const [copiedId, setCopiedId] = useState<string | null>(null)
  const t = useTranslations("common")

  if (!doc) {
    console.error("Doc not found in client component, should have been caught by server component.")
    notFound() // Or return a specific error component
  }

  const handleCopy = (text: string, id: string) => {
    navigator.clipboard
      .writeText(text)
      .then(() => {
        setCopiedId(id)
        setTimeout(() => setCopiedId(null), 2000) // Reset after 2 seconds
      })
      .catch((err) => {
        console.error("Failed to copy: ", err)
      })
  }

  // Define components object with proper typing for code
  const markdownComponents: Components = {
    // Add custom paragraph renderer
    p: ({ node, children }) => {
      // Check if the first child is our custom code block wrapper (a div)
      // or other block elements that shouldn't be nested in <p>
      if (node && node.children[0] && node.children[0].type === 'element' && node.children[0].tagName === 'div') {
        return <>{children}</> // Render children directly without <p>
      }
      // Otherwise, render the default paragraph
      return <p>{children}</p>
    },
    code: ({ node, inline, className, children, ...props }: CodeProps) => {
      // More robust check for inline code: check if the node's parent is NOT <pre>
      // The `inline` prop can sometimes be unreliable depending on plugins/markdown structure
      const isInline = !node?.position?.start.line === node?.position?.end.line || !className?.startsWith("language-");

      if (isInline) {
        // Render simple inline code
              return (
                <code className="bg-muted px-1.5 py-0.5 rounded text-sm font-medium text-foreground" {...props}>
                  {children}
                </code>
              )
            }

      // Handle block code (presumed if not inline)
            const match = /language-(\w+)/.exec(className || "")
            const language = match ? match[1] : ""
            const content = String(children)
            const codeId = `code-${content.length}-${content.slice(0, 20).replace(/\W/g, "")}`
            const isCopied = copiedId === codeId

            return (
              <div className="relative my-4">
                <pre className="rounded-md bg-slate-50 dark:bg-slate-800 p-4 overflow-x-auto">
                  <code className={`language-${language || "text"} text-foreground font-medium`} {...props}>
                    {children}
                  </code>
                </pre>
                {language && (
                  <div className="absolute top-2 left-3 text-xs font-medium text-slate-500 dark:text-slate-400">
                    {language === "r" || language === "R" ? "R" : language}
                  </div>
                )}
                <button
                  onClick={() => handleCopy(content, codeId)}
                  className="absolute top-2 right-2 rounded bg-slate-200/80 dark:bg-slate-700/80 p-1 text-xs font-medium text-slate-600 hover:bg-slate-300 dark:text-slate-300 dark:hover:bg-slate-600 transition-colors"
                  aria-label={isCopied ? "Copied" : "Copy code"}
                >
                  {isCopied ? (
                    <svg
                      xmlns="http://www.w3.org/2000/svg"
                      width="14"
                      height="14"
                      viewBox="0 0 24 24"
                      fill="none"
                      stroke="currentColor"
                      strokeWidth="2"
                      strokeLinecap="round"
                      strokeLinejoin="round"
                    >
                      <path d="M20 6L9 17l-5-5"></path>
                    </svg>
                  ) : (
                    <svg
                      xmlns="http://www.w3.org/2000/svg"
                      width="14"
                      height="14"
                      viewBox="0 0 24 24"
                      fill="none"
                      stroke="currentColor"
                      strokeWidth="2"
                      strokeLinecap="round"
                      strokeLinejoin="round"
                    >
                      <rect x="9" y="9" width="13" height="13" rx="2" ry="2"></rect>
                      <path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9a2 2 0 0 1 2 2v1"></path>
                    </svg>
                  )}
                </button>
              </div>
            )
          },
    h1: ({ node, ...props }) => <h1 id={slugify(props.children as string)} {...props} />,
    h2: ({ node, ...props }) => <h2 id={slugify(props.children as string)} {...props} />,
    h3: ({ node, ...props }) => <h3 id={slugify(props.children as string)} {...props} />,
    h4: ({ node, ...props }) => <h4 id={slugify(props.children as string)} {...props} />,
    h5: ({ node, ...props }) => <h5 id={slugify(props.children as string)} {...props} />,
    h6: ({ node, ...props }) => <h6 id={slugify(props.children as string)} {...props} />,
    img: ({ node, ...props }) => {
      // Check if alt text contains size information in format: alt|width=50%,height=50%
      const altTextParts = (props.alt || "").split('|');
      const altText = altTextParts[0];
      
      // Default styles with proper type definition
      let imgStyle: { [key: string]: string } = {};
      
      // Extract size information if present
      if (altTextParts.length > 1) {
        const sizeInfo = altTextParts[1];
        const sizeParts = sizeInfo.split(',');
        
        sizeParts.forEach(part => {
          const [key, value] = part.split('=');
          if (key && value) {
            imgStyle[key.trim()] = value.trim();
          }
        });
      }
      
      return (
            <img
              {...props}
          alt={altText}
          style={imgStyle}
              className="rounded-lg my-6 max-w-full"
              loading="lazy"
            />
      );
    },
  }

  return (
    <article className="prose prose-slate dark:prose-invert max-w-none">
      {doc.isFallback && (
        <div className="mb-4 p-3 bg-yellow-100 dark:bg-yellow-900/30 border border-yellow-300 dark:border-yellow-700 rounded-md text-sm text-yellow-800 dark:text-yellow-200">
          {t("translationNotAvailable")}
        </div>
      )}
      <h1>{doc.frontmatter.title}</h1>
      <ReactMarkdown
        rehypePlugins={[
          rehypeSlug,
          [
            rehypeAutolinkHeadings,
            {
              behavior: "wrap",
              properties: {
                className: ["anchor-link"],
              },
            },
          ],
        ]}
        components={markdownComponents} // Use the typed components object
      >
        {doc.content}
      </ReactMarkdown>
    </article>
  )
}
