"use client"

import { Bar, BarChart, CartesianGrid, Legend, ResponsiveContainer, Tooltip, XAxis, YAxis } from "recharts"

// Sample data for method comparison
const methodComparisonData = [
  {
    tissue: "Human Kidney",
    GPTcelltype: 0.78,
    Cassia: 0.85,
  },
  {
    tissue: "Human Lung",
    GPTcelltype: 0.82,
    Cassia: 0.87,
  },
  {
    tissue: "Human Large Intestine",
    GPTcelltype: 0.75,
    Cassia: 0.81,
  },
  {
    tissue: "Human Fetal Skin",
    GPTcelltype: 0.79,
    Cassia: 0.83,
  },
  {
    tissue: "Mouse Atlas",
    GPTcelltype: 0.73,
    Cassia: 0.8,
  },
]

export function MethodComparisonChart() {
  return (
    <div className="w-full h-[600px]">
      <ResponsiveContainer width="100%" height="100%">
        <BarChart
          data={methodComparisonData}
          margin={{
            top: 20,
            right: 30,
            left: 20,
            bottom: 70,
          }}
        >
          <CartesianGrid strokeDasharray="3 3" stroke="#f0f0f0" />
          <XAxis
            dataKey="tissue"
            angle={0}
            textAnchor="middle"
            height={60}
            tick={{ fill: "#4b5563" }}
            axisLine={{ stroke: "#e5e7eb" }}
            tickMargin={10}
          />
          <YAxis
            domain={[0, 1]}
            label={{
              value: "Accuracy Score",
              angle: -90,
              position: "insideLeft",
              style: { fill: "#4b5563", textAnchor: "middle" },
            }}
            tick={{ fill: "#4b5563" }}
            axisLine={{ stroke: "#e5e7eb" }}
            tickLine={{ stroke: "#e5e7eb" }}
          />
          <Tooltip
            formatter={(value, name) => [
              value.toFixed(2),
              name === "Cassia" ? "CASSIA (GPT-4o)" : "GPTcelltype (GPT-4o)",
            ]}
            contentStyle={{
              backgroundColor: "white",
              borderRadius: "0.5rem",
              boxShadow: "0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06)",
              border: "1px solid #e5e7eb",
            }}
            labelStyle={{ fontWeight: "bold", marginBottom: "0.5rem" }}
            itemStyle={{ padding: "2px 0" }}
            isAnimationActive={false}
          />
          <Legend
            verticalAlign="bottom"
            height={36}
            formatter={(value) => (
              <span style={{ color: "#4b5563", fontWeight: 500 }}>
                {value === "Cassia" ? "CASSIA (GPT-4o)" : "GPTcelltype (GPT-4o)"}
              </span>
            )}
          />
          <Bar dataKey="Cassia" name="Cassia" fill="#3b82f6" radius={[4, 4, 0, 0]} animationDuration={1500} />
          <Bar dataKey="GPTcelltype" name="GPTcelltype" fill="#10b981" radius={[4, 4, 0, 0]} animationDuration={1500} />
        </BarChart>
      </ResponsiveContainer>
    </div>
  )
}
