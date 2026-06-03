import json
from collections import Counter

path = '/home/runner/.claude/projects/-home-runner-work-shigella-shigella/2caf7287-97e8-4e7d-ac4a-b168a221f2d7.jsonl'
bash_cmds = Counter()
raw_cmds = []
mcp_tools = Counter()

with open(path) as f:
    for line in f:
        try:
            obj = json.loads(line)
        except Exception:
            continue
        msg = obj.get('message', {})
        if msg.get('role') != 'assistant':
            continue
        for block in msg.get('content', []):
            if not isinstance(block, dict) or block.get('type') != 'tool_use':
                continue
            name = block.get('name', '')
            inp = block.get('input', {})
            if name == 'Bash':
                cmd = inp.get('command', '').strip()
                raw_cmds.append(cmd)
                tokens = cmd.split()
                if tokens:
                    lead = tokens[0]
                    while '=' in lead and len(tokens) > 1:
                        tokens = tokens[1:]
                        lead = tokens[0]
                    sub = tokens[1] if len(tokens) > 1 else ''
                    key = (lead + ' ' + sub).strip()
                    bash_cmds[key] += 1
            elif name.startswith('mcp__'):
                mcp_tools[name] += 1

print('=== BASH (lead+sub) ===')
for k, v in bash_cmds.most_common(40):
    print(f'{v:3d}  {k}')

print()
print('=== RAW COMMANDS ===')
for cmd in raw_cmds:
    print(repr(cmd[:150]))

print()
print('=== MCP ===')
for k, v in mcp_tools.most_common(20):
    print(f'{v:3d}  {k}')
