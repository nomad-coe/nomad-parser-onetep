package eu.nomad_lab.parsers

import org.specs2.mutable.Specification

object OnetepParserSpec extends Specification {
  "OnetepParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(OnetepParser, "parsers/parser-onetep/test/examples/single_point.out", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(OnetepParser, "parsers/parser-onetep/test/examples/single_point.out", "json") must_== ParseResult.ParseSuccess
    }
  }
}

