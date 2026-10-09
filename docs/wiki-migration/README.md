# Non-destructive Wiki migration

Status: prepared for review; not published. The observed Wiki tip is
`715e5039c75e080814a12e957f5148c35cdf8bda`. Run the
[guarded preparation tool](../../tools/prepare_wiki.py) following the
[publishing instructions](../development/wiki.md).

| Page | Change | Canonical source |
| --- | --- | --- |
| Home | Replace navigation; preserve original as Historical Home | docs/index.md |
| Getting Started | Add concise entry | docs/getting-started/installation.md |
| Documentation | Add reference navigation | docs/api/index.md |
| Mathematical Conventions | Add convention entry | docs/architecture/conventions.md |
| Examples and Tutorials | Add learning paths | docs/examples/index.md and docs/tutorials/index.md |
| Development Roadmap | Link authoritative root roadmap | ROADMAP.md |
| Contributing | Add review workflow entry | docs/development/contributing.md |
| Historical and Legacy Notes | Explain historical status and compatibility | docs/api/legacy.md |
| Historical Home | Preserve original Home bytes and attribution | audit/wiki-snapshot/Home.md |
| _Sidebar | Add canonical navigation | Companion pages |
| Six original class/function pages | No changes | Existing Wiki / frozen snapshot |

[Home preview](pages/Home.md). `manifest.json` records original tip, planned paths
and file hashes. No source API manual is copied into the Wiki.

Historical Home retains the original trailing space and missing final newline to
preserve snapshot bytes; the whitespace diagnostic for that file is intentional.
