"use client"

import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select"

interface TissueSelectorProps {
  selectedTissue: string
  onTissueChange: (tissue: string) => void
}

export function TissueSelector({ selectedTissue, onTissueChange }: TissueSelectorProps) {
  return (
    <div className="w-full md:w-[250px]">
      <Select value={selectedTissue} onValueChange={onTissueChange}>
        <SelectTrigger className="tissue-selector bg-white border border-gray-200">
          <SelectValue placeholder="Select Tissue" />
        </SelectTrigger>
        <SelectContent>
          <SelectItem value="All Tissues">All Tissues</SelectItem>
          <SelectItem value="Human Kidney">Human Kidney</SelectItem>
          <SelectItem value="Human Lung">Human Lung</SelectItem>
          <SelectItem value="Human Large Intestine">Human Large Intestine</SelectItem>
          <SelectItem value="Human Fetal Skin">Human Fetal Skin</SelectItem>
          <SelectItem value="Mouse Atlas">Mouse Atlas</SelectItem>
        </SelectContent>
      </Select>
    </div>
  )
}
