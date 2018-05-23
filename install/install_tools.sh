#!/bin/bash
set -e

# RUN FROM REPOSITORY ROOT.

# Each directory in src/ is used to group tools by section.
# Tool man pages are defined in a flat directory.
# Every function has a Galaxy tool generated.
# Tools are defined in a configuration XML file which is deployed to Galaxy.

rm -rf tmp
mkdir -p tmp
cp install/stub_tool_conf.xml tmp/tool_conf.xml

find src -mindepth 1 -type f -name "*R" -print0 | while IFS= read -r -d $'\0' toolPath; do
  echo "Processing tool: $toolPath"

  groupName="$(dirname "$toolPath")"
  groupName=${groupName##*/}

  toolFile=$(basename "$toolPath")
  toolName="${toolFile%.*}"
  echo " - installing tool '$toolName' in group: $groupName"

  Rscript install/install_tool.R "$toolName" "$groupName"

  # Can optionally patch the tool XML
  # e.g. to support command lines which RGalaxy cannot yet generate.
  toolDiff=$(dirname "$toolPath")/${toolName}.xml.diff
  if [ -f "$toolDiff" ]; then
    echo " - Found tool XML diff at: $toolDiff"
    patch tmp/tools/${toolName}/${toolName}.xml "$toolDiff"
  fi

done

# Deploy the generated tools.
scp tmp/tool_conf.xml galaxy@localhost:installations/galaxy-dist/config/galaxy_tool_conf.xml
scp -r tmp/tools/* galaxy@localhost:installations/rgalaxy-tools/

./install/install_alt_tools.sh