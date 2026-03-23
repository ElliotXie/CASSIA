import type { Metadata } from 'next';
import AnnotationBoostClient from './AnnotationBoostClient';

export const metadata: Metadata = {
  title: 'Annotation Boost Agent - CASSIA',
  description: 'Iterative marker analysis for deeper cell type annotation with AI.',
};

export default function AnnotationBoostPage() {
  return <AnnotationBoostClient />;
}
