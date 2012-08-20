#!/usr/bin/env bash

LOC="lma27ad.f"
TMP=$(mktemp)

cd $(dirname "$0")
cat "$1" > "$LOC"

BEGIN=$(egrep -in 'SUBROUTINE\ +MA27OD' "$LOC" | cut -d: -f1)

for LINE in $(egrep -in 'INFO\(1\)\ +=\ +-5' "$LOC" | cut -d: -f1); do
  if [ $LINE -gt $BEGIN ]; then
    head -n  $((LINE+1)) "$LOC"   > "$TMP"
    echo '      INFO(16) = APOS' >> "$TMP"
    tail -n +$((LINE+2)) "$LOC"  >> "$TMP"
    cat "$TMP"                    > "$LOC"
    break
  fi
done

for LINE in $(egrep -in 'INFO\(1\)\ +=\ +-6' "$LOC" | cut -d: -f1); do
  if [ $LINE -gt $BEGIN ]; then
    head -n  $((LINE+1)) "$LOC"   > "$TMP"
    echo '      INFO(16) = APOS' >> "$TMP"
    tail -n +$((LINE+2)) "$LOC"  >> "$TMP"
    cat "$TMP"                    > "$LOC"
    break
  fi
done

for LINE in $(egrep -in 'INFO\(1\)\ +=\ +3' "$LOC" | cut -d: -f1); do
  if [ $LINE -gt $BEGIN ]; then
    head -n  $((LINE+1)) "$LOC"   > "$TMP"
    echo '      INFO(16) = APOS' >> "$TMP"
    tail -n +$((LINE+2)) "$LOC"  >> "$TMP"
    cat "$TMP"                    > "$LOC"
    break
  fi
done

for LINE in $(egrep -in 'INFO\(1\)\ +=\ +2' "$LOC" | cut -d: -f1); do
  if [ $LINE -gt $BEGIN ]; then
    head -n $LINE "$LOC"                > "$TMP"
    echo '            INFO(16) = APOS' >> "$TMP"
    tail -n +$((LINE+1)) "$LOC"        >> "$TMP"
    cat "$TMP"                          > "$LOC"
    break
  fi
done

rm -f "$TMP"
exit 0
