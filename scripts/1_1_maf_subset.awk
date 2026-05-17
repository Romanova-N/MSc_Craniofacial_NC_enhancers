NR==FNR { keep[$1]=1; next }

function flush() {
  if (s_count >= 1) for (i=1; i<=n; i++) print buf[i]
  n=0; s_count=0
}

BEGIN { n=0; s_count=0; keep["hg38"]=1 }

/^##maf/ { print; next }
/^a[ \t]/ {
  if (n>0) flush()
  buf[++n]=$0
  next
}
/^s[ \t]/ {
  split($2, a, "."); genome=a[1]
  if (keep[genome]) { buf[++n]=$0; s_count++ }
  next
}
/^[ \t]*$/ {
  if (n>0) { buf[++n]=""; flush() }
  next
}

END { if (n>0) flush() }