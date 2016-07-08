package eu.nomad_lab.parsers

import org.specs2.mutable.Specification

object OnetepParserSpec extends Specification {
  "OnetepParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(OnetepParser, "parsers/onetep/test/examples/test2/single_point_2.out", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(OnetepParser, "parsers/onetep/test/examples/test2/single_point_2.out", "json") must_== ParseResult.ParseSuccess
    }
  }
}

