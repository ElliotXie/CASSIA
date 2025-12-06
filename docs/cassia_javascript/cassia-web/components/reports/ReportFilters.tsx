'use client'

import React, { useState, useMemo } from 'react'
import { Search, X, Filter } from 'lucide-react'
import { Input } from '@/components/ui/input'
import { Button } from '@/components/ui/button'
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from '@/components/ui/select'
import { cn } from '@/lib/utils'

export interface FilterOption {
  value: string
  label: string
}

export interface FilterConfig {
  key: string
  label: string
  options: FilterOption[]
  placeholder?: string
}

export interface ReportFiltersProps {
  searchPlaceholder?: string
  searchValue: string
  onSearchChange: (value: string) => void
  filters?: FilterConfig[]
  filterValues?: Record<string, string>
  onFilterChange?: (key: string, value: string) => void
  onClearAll?: () => void
  className?: string
  showStats?: boolean
  totalCount?: number
  visibleCount?: number
}

/**
 * Reusable search and filter controls for reports
 * Supports text search and multiple dropdown filters
 */
export function ReportFilters({
  searchPlaceholder = 'Search clusters, cell types, markers...',
  searchValue,
  onSearchChange,
  filters = [],
  filterValues = {},
  onFilterChange,
  onClearAll,
  className,
  showStats = true,
  totalCount = 0,
  visibleCount = 0
}: ReportFiltersProps) {
  const hasActiveFilters = useMemo(() => {
    return searchValue.length > 0 || Object.values(filterValues).some(v => v && v !== '')
  }, [searchValue, filterValues])

  return (
    <div className={cn(
      "bg-white/80 backdrop-blur-sm border border-teal-100 rounded-xl p-4 shadow-sm",
      className
    )}>
      <div className="flex flex-wrap gap-3 items-center">
        {/* Search Input */}
        <div className="relative flex-1 min-w-[250px]">
          <Search className="absolute left-3 top-1/2 -translate-y-1/2 h-4 w-4 text-gray-400" />
          <Input
            type="text"
            placeholder={searchPlaceholder}
            value={searchValue}
            onChange={(e) => onSearchChange(e.target.value)}
            className="pl-10 bg-white border-gray-200 focus:border-teal-500 focus:ring-teal-500/20"
          />
        </div>

        {/* Filter Dropdowns */}
        {filters.map((filter) => (
          <Select
            key={filter.key}
            value={filterValues[filter.key] || ''}
            onValueChange={(value) => onFilterChange?.(filter.key, value)}
          >
            <SelectTrigger className="w-[140px] bg-white border-gray-200">
              <SelectValue placeholder={filter.placeholder || filter.label} />
            </SelectTrigger>
            <SelectContent>
              <SelectItem value="">All {filter.label}</SelectItem>
              {filter.options.map((option) => (
                <SelectItem key={option.value} value={option.value}>
                  {option.label}
                </SelectItem>
              ))}
            </SelectContent>
          </Select>
        ))}

        {/* Clear Filters Button */}
        {hasActiveFilters && onClearAll && (
          <Button
            variant="outline"
            size="sm"
            onClick={onClearAll}
            className="border-gray-200 hover:border-teal-500 hover:text-teal-600"
          >
            <X className="h-4 w-4 mr-1" />
            Clear
          </Button>
        )}
      </div>

      {/* Stats Display */}
      {showStats && (
        <div className="flex items-center gap-2 mt-3 text-sm text-gray-500">
          <Filter className="h-3.5 w-3.5" />
          <span>
            Showing <span className="font-semibold text-teal-700">{visibleCount}</span>
            {totalCount !== visibleCount && (
              <> of <span className="font-semibold">{totalCount}</span></>
            )}
            {totalCount === 1 ? ' cluster' : ' clusters'}
          </span>
        </div>
      )}
    </div>
  )
}

/**
 * Hook to manage filter state for reports
 */
export function useReportFilters<T extends Record<string, any>>(
  data: T[],
  searchFields: (keyof T)[],
  filterFields?: (keyof T)[]
) {
  const [searchValue, setSearchValue] = useState('')
  const [filterValues, setFilterValues] = useState<Record<string, string>>({})

  // Extract unique filter options from data
  const filterOptions = useMemo(() => {
    const options: Record<string, FilterOption[]> = {}

    filterFields?.forEach((field) => {
      const uniqueValues = [...new Set(
        data
          .map(item => String(item[field] || ''))
          .filter(v => v && v !== 'N/A' && v !== 'undefined' && v !== 'null')
      )].sort()

      options[String(field)] = uniqueValues.map(v => ({ value: v, label: v }))
    })

    return options
  }, [data, filterFields])

  // Filter data based on search and filters
  const filteredData = useMemo(() => {
    return data.filter(item => {
      // Search filter
      if (searchValue) {
        const searchLower = searchValue.toLowerCase()
        const matchesSearch = searchFields.some(field => {
          const value = String(item[field] || '').toLowerCase()
          return value.includes(searchLower)
        })
        if (!matchesSearch) return false
      }

      // Dropdown filters
      for (const [key, value] of Object.entries(filterValues)) {
        if (value && value !== '') {
          if (String(item[key as keyof T]) !== value) {
            return false
          }
        }
      }

      return true
    })
  }, [data, searchValue, filterValues, searchFields])

  const handleFilterChange = (key: string, value: string) => {
    setFilterValues(prev => ({ ...prev, [key]: value }))
  }

  const handleClearAll = () => {
    setSearchValue('')
    setFilterValues({})
  }

  return {
    searchValue,
    setSearchValue,
    filterValues,
    setFilterValues,
    filterOptions,
    filteredData,
    handleFilterChange,
    handleClearAll,
    totalCount: data.length,
    visibleCount: filteredData.length
  }
}

export default ReportFilters
