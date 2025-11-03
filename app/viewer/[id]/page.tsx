import { redirect } from "next/navigation";

interface ViewerIdPageProps {
  params: Promise<{ id: string }>;
}

export default async function ViewerIdRedirect({ params }: ViewerIdPageProps) {
  const { id } = await params;
  const search = new URLSearchParams({ dataset: id });

  redirect(`/?${search.toString()}`);
}
