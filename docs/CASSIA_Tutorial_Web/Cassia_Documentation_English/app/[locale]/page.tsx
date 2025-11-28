import Link from "next/link"
import { Button } from "@/components/ui/button"
import { Dna, Microscope, FileSearch, CheckCircle } from "lucide-react"
import { useTranslations } from "next-intl"
import { LanguageSwitcher } from "@/components/language-switcher"

export default function Home() {
  const t = useTranslations("home")

  return (
    <div className="flex min-h-screen flex-col">
      {/* Language switcher header */}
      <div className="flex justify-end items-center gap-3 p-4">
        <span className="text-sm text-muted-foreground">{t("language")}</span>
        <LanguageSwitcher />
      </div>
      <div className="flex flex-col items-center justify-center px-4 py-16 text-center md:py-28">
        <div className="flex h-20 w-20 items-center justify-center rounded-2xl bg-primary text-primary-foreground mb-6">
          <Dna className="h-12 w-12" />
        </div>
        <h1 className="text-4xl font-extrabold tracking-tight sm:text-5xl md:text-6xl lg:text-7xl">{t("title")}</h1>
        <p className="mt-4 text-xl text-muted-foreground max-w-2xl">
          {t("subtitle")}
        </p>
        <div className="mt-8 flex flex-wrap justify-center gap-4">
          <Button asChild variant="outline" size="lg" className="h-12 px-8">
            <Link href="https://github.com/ElliotXie/CASSIA" target="_blank" rel="noopener noreferrer">
              {t("viewOnGithub")}
            </Link>
          </Button>
          <Button asChild variant="outline" size="lg" className="h-12 px-8">
            <Link
              href="https://www.biorxiv.org/content/10.1101/2024.12.04.626476v2"
              target="_blank"
              rel="noopener noreferrer"
            >
              {t("readThePaper")}
            </Link>
          </Button>
        </div>
      </div>

      <div className="container mx-auto px-4 py-12 md:py-24">
        <div className="grid gap-8 md:grid-cols-3">
          <div className="rounded-lg border bg-card p-6 shadow-sm transition-all hover:shadow-md">
            <div className="mb-4 flex h-10 w-10 items-center justify-center rounded-full bg-primary/10 text-primary">
              <Microscope className="h-5 w-5" />
            </div>
            <h3 className="text-lg font-semibold">{t("interpretable")}</h3>
            <p className="mt-2 text-muted-foreground">
              {t("interpretableDesc")}
            </p>
          </div>
          <div className="rounded-lg border bg-card p-6 shadow-sm transition-all hover:shadow-md">
            <div className="mb-4 flex h-10 w-10 items-center justify-center rounded-full bg-primary/10 text-primary">
              <FileSearch className="h-5 w-5" />
            </div>
            <h3 className="text-lg font-semibold">{t("referenceFree")}</h3>
            <p className="mt-2 text-muted-foreground">
              {t("referenceFreeDesc")}
            </p>
          </div>
          <div className="rounded-lg border bg-card p-6 shadow-sm transition-all hover:shadow-md">
            <div className="mb-4 flex h-10 w-10 items-center justify-center rounded-full bg-primary/10 text-primary">
              <CheckCircle className="h-5 w-5" />
            </div>
            <h3 className="text-lg font-semibold">{t("accurate")}</h3>
            <p className="mt-2 text-muted-foreground">
              {t("accurateDesc")}
            </p>
          </div>
        </div>
      </div>
    </div>
  )
}
