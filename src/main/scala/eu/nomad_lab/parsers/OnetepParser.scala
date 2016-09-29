package eu.nomad_lab.parsers

import eu.{ nomad_lab => lab }
import eu.nomad_lab.DefaultPythonInterpreter
import org.{ json4s => jn }
import scala.collection.breakOut

object OnetepParser extends SimpleExternalParserGenerator(
  name = "OnetepParser",
  parserInfo = jn.JObject(
    ("name" -> jn.JString("OnetepParser")) ::
      ("parserId" -> jn.JString("OnetepParser" + lab.OnetepVersionInfo.version)) ::
      ("versionInfo" -> jn.JObject(
        ("nomadCoreVersion" -> jn.JObject(lab.NomadCoreVersionInfo.toMap.map {
          case (k, v) => k -> jn.JString(v.toString)
        }(breakOut): List[(String, jn.JString)])) ::
          (lab.OnetepVersionInfo.toMap.map {
            case (key, value) =>
              (key -> jn.JString(value.toString))
          }(breakOut): List[(String, jn.JString)])
      )) :: Nil
  ),
  mainFileTypes = Seq("text/.*"),
  mainFileRe = """\s*\|\s*Linear-Scaling Ab Initio Total Energy Program\s*\|\s*""".r,
  cmd = Seq(DefaultPythonInterpreter.pythonExe(), "${envDir}/parsers/onetep/parser/parser-onetep/OnetepParser.py",
    "--uri", "${mainFileUri}", "${mainFilePath}"),
  resList = Seq(
    "parser-onetep/OnetepParser.py",
    "parser-onetep/OnetepBandParser.py",
    "parser-onetep/OnetepCellParser.py",
    "parser-onetep/OnetepCommon.py",
    "parser-onetep/OnetepMDParser.py",
    "parser-onetep/OnetepTSParser.py",
    "parser-onetep/setup_paths.py",
    "nomad_meta_info/public.nomadmetainfo.json",
    "nomad_meta_info/common.nomadmetainfo.json",
    "nomad_meta_info/meta_types.nomadmetainfo.json",
    "nomad_meta_info/onetep.nomadmetainfo.json"
  ) ++ DefaultPythonInterpreter.commonFiles(),
  dirMap = Map(
    "parser-onetep" -> "parsers/onetep/parser/parser-onetep",
    "nomad_meta_info" -> "nomad-meta-info/meta_info/nomad_meta_info"
  ) ++ DefaultPythonInterpreter.commonDirMapping()
)
