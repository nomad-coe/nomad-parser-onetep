/*
 * Copyright 2016-2018 Martina Stella, Fawzi Mohamed
 * 
 *   Licensed under the Apache License, Version 2.0 (the "License");
 *   you may not use this file except in compliance with the License.
 *   You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 *   Unless required by applicable law or agreed to in writing, software
 *   distributed under the License is distributed on an "AS IS" BASIS,
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *   See the License for the specific language governing permissions and
 *   limitations under the License.
 */

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

object OnetepParserSpec_1 extends Specification {
  "OnetepParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(OnetepParser, "parsers/onetep/test/examples/test08/ethene_relax.out", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(OnetepParser, "parsers/onetep/test/examples/test08/ethene_relax.out", "json") must_== ParseResult.ParseSuccess
    }
  }
}

object OnetepParserSpec_2 extends Specification {
  "OnetepParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(OnetepParser, "parsers/onetep/test/examples/test49/test_49.out", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(OnetepParser, "parsers/onetep/test/examples/test49/test_49.out", "json") must_== ParseResult.ParseSuccess
    }
  }
}
object OnetepParserSpec_3 extends Specification {
  "OnetepParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(OnetepParser, "parsers/onetep/test/examples/test20/test20.out", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(OnetepParser, "parsers/onetep/test/examples/test20/test20.out", "json") must_== ParseResult.ParseSuccess
    }
  }
}
object OnetepParserSpec_4 extends Specification {
  "OnetepParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(OnetepParser, "parsers/onetep/test/examples/ethene/ethene_tddft.out", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(OnetepParser, "parsers/onetep/test/examples/ethene/ethene_tddft.out", "json") must_== ParseResult.ParseSuccess
    }
  }
}
object OnetepParserSpec_5 extends Specification {
  "OnetepParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(OnetepParser, "parsers/onetep/test/examples/TS/reac.out", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(OnetepParser, "parsers/onetep/test/examples/TS/reac.out", "json") must_== ParseResult.ParseSuccess
    }
  }
}
object OnetepParserSpec_6 extends Specification {
  "OnetepParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(OnetepParser, "parsers/onetep/test/examples/benz_dim/benzene_dimer_vdw.out", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(OnetepParser, "parsers/onetep/test/examples/benz_dim/benzene_dimer_vdw.out", "json") must_== ParseResult.ParseSuccess
    }
  }
}
object OnetepParserSpec_7 extends Specification {
  "OnetepParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(OnetepParser, "parsers/onetep/test/examples/fluor/12-difluoroethane.out", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(OnetepParser, "parsers/onetep/test/examples/fluor/12-difluoroethane.out", "json") must_== ParseResult.ParseSuccess
    }
  }
}