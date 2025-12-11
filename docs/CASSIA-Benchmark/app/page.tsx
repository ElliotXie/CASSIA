import { redirect } from "next/navigation"
import Link from "next/link"
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card"
import { MethodComparisonChart } from "@/components/method-comparison-chart"
import { Button } from "@/components/ui/button"
import { ChevronRight, BarChart2, FileText, Database } from "lucide-react"

export default function Home() {
  // Redirect to the CASSIA benchmark page
  redirect("/methods/cassia")
}
