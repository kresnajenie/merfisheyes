import { redirect } from "next/navigation";

interface ViewerIdPageProps {
  params: { id: string };
}

export default function ViewerIdRedirect({ params }: ViewerIdPageProps) {
  const search = new URLSearchParams({ dataset: params.id });

  redirect(`/?${search.toString()}`);
}
