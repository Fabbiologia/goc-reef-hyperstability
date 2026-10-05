#!/usr/bin/env python3
"""Post-process a pandoc-generated .docx: A4 page setup, continuous line
numbers, and a right-aligned page-number footer (no running head). Usage: python docx_linenums_pagenums.py file.docx"""
import sys, re, zipfile, shutil, os

RNS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
FOOTER_CT = "application/vnd.openxmlformats-officedocument.wordprocessingml.footer+xml"

NEW_SECTPR = (
    '<w:sectPr>'
    '<w:footerReference w:type="default" r:id="rId900"/>'
    '<w:footnotePr><w:numRestart w:val="eachSect"/></w:footnotePr>'
    '<w:pgSz w:w="11906" w:h="16838"/>'
    '<w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440" '
    'w:header="708" w:footer="708" w:gutter="0"/>'
    '<w:lnNumType w:countBy="1" w:restart="continuous" w:distance="284"/>'
    '<w:cols w:space="708"/>'
    '</w:sectPr>'
)

FOOTER_XML = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
    '<w:ftr xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" '
    'xmlns:r="' + RNS + '">'
    '<w:p><w:pPr><w:jc w:val="right"/></w:pPr>'
    '<w:r><w:fldChar w:fldCharType="begin"/></w:r>'
    '<w:r><w:instrText xml:space="preserve"> PAGE </w:instrText></w:r>'
    '<w:r><w:fldChar w:fldCharType="end"/></w:r>'
    '</w:p></w:ftr>'
)

def patch(path):
    tmp = path + ".tmp"
    zin = zipfile.ZipFile(path)
    names = zin.namelist()
    doc = zin.read("word/document.xml").decode("utf-8")
    rels = zin.read("word/_rels/document.xml.rels").decode("utf-8")
    ct = zin.read("[Content_Types].xml").decode("utf-8")

    # 1) replace the body sectPr (the last one) with a full sectPr
    sectprs = list(re.finditer(r'<w:sectPr\b.*?</w:sectPr>', doc, re.S))
    if not sectprs:
        raise SystemExit("no sectPr found")
    last = sectprs[-1]
    doc = doc[:last.start()] + NEW_SECTPR + doc[last.end():]
    # ensure xmlns:r is declared on the <w:document> element (pandoc usually does)
    root = re.search(r'<w:document\b[^>]*>', doc).group(0)
    if "xmlns:r=" not in root:
        doc = doc.replace(root, root[:-1] + ' xmlns:r="%s">' % RNS, 1)

    # 2) footer relationship
    if 'Id="rId900"' not in rels:
        rels = rels.replace("</Relationships>",
            '<Relationship Id="rId900" Type="%s/footer" Target="footer1.xml"/></Relationships>' % RNS)

    # 3) content-type override for the footer part
    if "footer1.xml" not in ct:
        ct = ct.replace("</Types>",
            '<Override PartName="/word/footer1.xml" ContentType="%s"/></Types>' % FOOTER_CT)

    # write everything back
    zout = zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED)
    for n in names:
        if n == "word/document.xml":
            zout.writestr(n, doc)
        elif n == "word/_rels/document.xml.rels":
            zout.writestr(n, rels)
        elif n == "[Content_Types].xml":
            zout.writestr(n, ct)
        else:
            zout.writestr(n, zin.read(n))
    if "word/footer1.xml" not in names:
        zout.writestr("word/footer1.xml", FOOTER_XML)
    zin.close(); zout.close()
    shutil.move(tmp, path)
    print("patched %s: line numbers + page-number footer + A4 page setup" % os.path.basename(path))

if __name__ == "__main__":
    patch(sys.argv[1])
