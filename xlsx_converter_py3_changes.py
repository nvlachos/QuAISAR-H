--- xlsx_converter.py	(original)
+++ xlsx_converter.py	(refactored)
@@ -49,7 +49,7 @@
 elif not sys.stdin.isatty():
   z = zipfile.ZipFile(sys.stdin)
 else:
-  print __doc__.decode().strip()
+  print(__doc__.decode().strip())
   sys.exit(1)

 n=lambda x: "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}%s" % x
@@ -60,13 +60,13 @@

 def sheet_report():
   global sheet_xs
-  print>>sys.stderr, "Sheets in this file:"
+  print("Sheets in this file:", file=sys.stderr)
   for i,x in enumerate(sheet_xs):
-    print>>sys.stderr, "%3d: %s" % (i+1, x.get('name'))
+    print("%3d: %s" % (i+1, x.get('name')), file=sys.stderr)
   sys.exit(1)

 def sheet_error(msg):
-  print>>sys.stderr, msg
+  print(msg, file=sys.stderr)
   sheet_report()

 if not args and len(sheet_filenames) > 1:
@@ -134,7 +134,7 @@
   global warning_count
   warning_count += 1
   if warning_count > warning_max: return
-  print>>sys.stderr, "WARNING: %s" % s
+  print("WARNING: %s" % s, file=sys.stderr)

 def cell_text_clean(text):
   s = text.encode("utf-8")
@@ -154,7 +154,7 @@
   for c,j in zip(cells_elts,inds):
     cells[j] = c
   #print( *(cell2text( c ).encode("utf-8").replace("\t"," ") for c in cells), sep="\t")
-  print myjoin((cell_text_clean(cell2text( c )) for c in cells), sep="\t")
+  print(myjoin((cell_text_clean(cell2text( c )) for c in cells), sep="\t"))

 if warning_count > warning_max:
-  print>>sys.stderr, "%d total warnings, %d hidden" % (warning_count, warning_count-warning_max)
+  print("%d total warnings, %d hidden" % (warning_count, warning_count-warning_max), file=sys.stderr)
