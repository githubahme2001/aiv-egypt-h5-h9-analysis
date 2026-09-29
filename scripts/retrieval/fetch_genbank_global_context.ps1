<#
  Global contextual H5 haemagglutinin sequences for the Egyptian 2.3.4.4b phylogeny.

  Purpose: the Egypt-only tree cannot distinguish "the post-2021 Egyptian H5N1 descends from the
  resident Egyptian 2.3.4.4b population" from "it is a separate wild-bird introduction of the same
  clade" (Naguib et al., Pathogens 2023;12:90). Adding contextual taxa from the countries along the
  same flyway, and from the 2016-2019 H5N8 panzootic, makes that testable.

  Output lands next to this script. Same cut-off as the main dataset so the two are comparable.

  Run:
      cd "scripts\retrieval"
      .\fetch_genbank_global_context.ps1
#>

[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12

$outDir  = if ($PSScriptRoot) { $PSScriptRoot } else { (Get-Location).Path }
$base    = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
$email   = "your.email@example.org"
$extra   = "&tool=EgyptAIV_context&email=$([uri]::EscapeDataString($email))"
$CUTOFF  = "2026/09/28"
$utf8    = New-Object System.Text.UTF8Encoding($false)
$logPath = Join-Path $outDir "NCBI_context_log.txt"

Write-Host ""
Write-Host "  Global contextual H5 pull  ->  $outDir" -ForegroundColor Cyan
Write-Host ""

function Write-Log($m) { Write-Host $m; Add-Content -Path $logPath -Value $m -Encoding UTF8 }
function Get-NodeText($node, $name) {
    if (-not $node) { return "" }
    $n = $node.SelectSingleNode($name)
    if ($n) { return $n.InnerText } else { return "" }
}
function Invoke-Eutil($url) {
    for ($try = 1; $try -le 4; $try++) {
        try { return Invoke-WebRequest -Uri $url -UseBasicParsing -TimeoutSec 120 }
        catch {
            if ($try -eq 4) { throw }
            Write-Log "    retry $try after error: $($_.Exception.Message)"
            Start-Sleep -Seconds (3 * $try)
        }
    }
}

# Countries chosen for the Black Sea / Mediterranean and East African flyways, plus the
# 2016-2019 H5N8 panzootic range. Egypt is deliberately EXCLUDED - it is already in the main set.
$countries = @(
    "Nigeria","Niger","Ghana","Cameroon","Burkina Faso","South Africa","Lesotho",
    "Israel","Iraq","Saudi Arabia","Kuwait","Iran","Turkey",
    "Netherlands","Germany","United Kingdom","France","Italy","Poland","Hungary",
    "Russia","Kazakhstan","Ukraine","Romania","Bulgaria","Greece"
)

$geo = ($countries | ForEach-Object { "`"$_`"[All Fields]" }) -join " OR "
$h5  = @("H5N1","H5N2","H5N5","H5N6","H5N8","H5Nx","H5") | ForEach-Object { "`"$_`"[All Fields]" }
$h5q = $h5 -join " OR "
$gen = @("hemagglutinin","haemagglutinin","HA") | ForEach-Object { "`"$_`"[All Fields]" }
$gnq = $gen -join " OR "
$cut = "(`"2015/01/01`"[PDAT] : `"$CUTOFF`"[PDAT])"

$query = "(`"Influenza A virus`"[Organism]) AND ($geo) AND ($h5q) AND (($gnq) OR (`"segment 4`"[All Fields])) AND $cut NOT (Egypt[All Fields])"

Write-Log "=== Global contextual H5 HA, 26 countries, 2015-$CUTOFF ==="
Write-Log "  query: $query"

$u = "$base/esearch.fcgi?db=nuccore&term=$([uri]::EscapeDataString($query))&usehistory=y&retmax=0$extra"
$doc = New-Object System.Xml.XmlDocument
$doc.LoadXml((Invoke-Eutil $u).Content)
$root   = $doc.DocumentElement
$count  = [int](Get-NodeText $root "Count")
$webenv = Get-NodeText $root "WebEnv"
$qkey   = Get-NodeText $root "QueryKey"
if (-not $webenv) { throw "esearch returned no history handle" }
Write-Log "  $count records found"

if ($count -gt 20000) {
    Write-Log "  NOTE: $count is large. Consider narrowing the country list before proceeding."
}

$fastaOut = Join-Path $outDir "NCBI_CONTEXT_H5_HA.fasta"
$metaOut  = Join-Path $outDir "NCBI_CONTEXT_H5_HA_metadata.tsv"

# ---- sequences ----
$sb = New-Object System.Text.StringBuilder
for ($start = 0; $start -lt $count; $start += 200) {
    $u = "$base/efetch.fcgi?db=nuccore&query_key=$qkey&WebEnv=$webenv&retstart=$start&retmax=200&rettype=fasta&retmode=text$extra"
    [void]$sb.Append((Invoke-Eutil $u).Content)
    Write-Log ("    sequences {0}/{1}" -f ([Math]::Min($start + 200, $count)), $count)
    Start-Sleep -Milliseconds 400
}
[System.IO.File]::WriteAllText($fastaOut, $sb.ToString(), $utf8)
Write-Log "  wrote sequences -> $fastaOut"

# ---- metadata ----
$header = "accession`tlength`tcountry`thost`tcollection_date`tserotype`tsegment`tstrain`ttitle"
$rows = New-Object System.Collections.Generic.List[string]
$rows.Add($header)
for ($start = 0; $start -lt $count; $start += 200) {
    $u = "$base/esummary.fcgi?db=nuccore&query_key=$qkey&WebEnv=$webenv&retstart=$start&retmax=200&version=2.0$extra"
    $sdoc = New-Object System.Xml.XmlDocument
    $sdoc.LoadXml((Invoke-Eutil $u).Content)
    foreach ($ds in $sdoc.GetElementsByTagName("DocumentSummary")) {
        $acc   = Get-NodeText $ds "AccessionVersion"
        $len   = Get-NodeText $ds "Slen"
        $title = Get-NodeText $ds "Title"
        $sub   = Get-NodeText $ds "SubType"
        $val   = Get-NodeText $ds "SubName"
        $country=""; $host=""; $cdate=""; $sero=""; $seg=""; $strain=""
        if ($sub -and $val) {
            $ks = $sub -split "\|"; $vs = $val -split "\|"
            for ($i = 0; $i -lt [Math]::Min($ks.Count, $vs.Count); $i++) {
                switch ($ks[$i]) {
                    "country"         { $country = $vs[$i] }
                    "geo_loc_name"    { if (-not $country) { $country = $vs[$i] } }
                    "host"            { $host    = $vs[$i] }
                    "collection_date" { $cdate   = $vs[$i] }
                    "serotype"        { $sero    = $vs[$i] }
                    "segment"         { $seg     = $vs[$i] }
                    "strain"          { $strain  = $vs[$i] }
                    "isolate"         { if (-not $strain) { $strain = $vs[$i] } }
                }
            }
        }
        $rows.Add(($acc, $len, $country, $host, $cdate, $sero, $seg, $strain, $title) -join "`t")
    }
    Write-Log ("    metadata  {0}/{1}" -f ([Math]::Min($start + 200, $count)), $count)
    Start-Sleep -Milliseconds 400
}
[System.IO.File]::WriteAllLines($metaOut, $rows, $utf8)
Write-Log "  wrote metadata  -> $metaOut"

Write-Log ""
Write-Log "=== done: $count contextual records ==="
Write-Host ""
Write-Host "  Finished. Tell Claude the context files are ready." -ForegroundColor Green
Write-Host ""
