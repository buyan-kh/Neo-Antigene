/**
 * Empty states always say why they are empty.
 *
 * An empty panel that looks like a loading state, or worse one padded with
 * illustrative numbers, is the fastest way to lose a technical audience.
 */
export function EmptyState({
  title,
  body,
  action,
}: {
  title: string;
  body: string;
  action?: React.ReactNode;
}) {
  return (
    <div className="mx-auto flex max-w-lg flex-col items-center gap-3 rounded-lg border border-dashed px-6 py-16 text-center">
      <h2 className="text-sm font-medium">{title}</h2>
      <p className="text-sm leading-relaxed text-muted-foreground">{body}</p>
      {action && <div className="mt-2">{action}</div>}
    </div>
  );
}
