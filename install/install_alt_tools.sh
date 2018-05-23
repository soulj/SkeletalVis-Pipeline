#!/bin/bash
set -e

# RUN FROM REPOSITORY ROOT
# Run automatically by install_tools.sh during installation

echo "Installing Kallisto tools."

for xmlfile in alt-src/kallisto/*.xml; do
  FILENAME=$(basename "$xmlfile" .xml)
  if [ "$FILENAME" != "kallisto_tool_conf" ]; then
    rsync -az $xmlfile galaxy@localhost:installations/rgalaxy-tools/"$FILENAME"/ 
  fi
done

for rfile in alt-src/kallisto/*.R; do
  FILENAME=$(basename "$rfile" .R)
  rsync -az $rfile galaxy@localhost:installations/rgalaxy-tools/"$FILENAME"/
done  

echo "Copying kallisto tool panel config."

scp alt-src/kallisto/kallisto_tool_conf.xml galaxy@localhost:installations/galaxy-dist/config/kallisto_tool_conf.xml

echo "Kallisto tool installation complete."

echo "Installing RNASeqExpressionData tool."

TOOL="getRNASeqExpressionData"
TOOLDIR="alt-src/getRNASeqExpressionData/"
rsync -az $TOOLDIR galaxy@localhost:installations/rgalaxy-tools/"$TOOL"/

# Add to tool panel config
xmlstarlet ed -s "/toolbox/section[@name='SYBIL Systems Biology']" -t elem -n tool -v "" tmp/tool_conf.xml > tmp/alt_tool_conf.xml
xmlstarlet ed --inplace -i "/toolbox/section[@name='SYBIL Systems Biology']/tool[not(@file)]" -t attr -n file -v "getRNASeqExpressionData/getRNASeqExpressionData.xml" tmp/alt_tool_conf.xml
scp tmp/alt_tool_conf.xml galaxy@localhost:installations/galaxy-dist/config/galaxy_tool_conf.xml
echo "RNASeqExpressionData tool install complete."

echo "Installing QCTools"
TOOLDIR="alt-src/QCTools/"
for directory in $TOOLDIR*/; do
  DIRNAME=$(basename "$directory")
  rsync -az $directory galaxy@localhost:installations/rgalaxy-tools/"$DIRNAME"/
done
scp ${TOOLDIR}qctools_tool_conf.xml galaxy@localhost:installations/galaxy-dist/config/qctools_tools_conf.xml

echo "QCTools installation complete."