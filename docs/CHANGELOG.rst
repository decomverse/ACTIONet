*********
Changelog
*********

0.4.13
======
:Date: May 18, 2021

* `220 <https://github.com/insitro/redun/pull/220>`_ - DE-2637 fix hashing of task_options_update
* `204 <https://github.com/insitro/redun/pull/204>`_ - DE-2619 Use O(1) queries to speedup record serialization
* `218 <https://github.com/insitro/redun/pull/218>`_ - DE-2635 Show unknown CallNodes for unfinished jobs
* `217 <https://github.com/insitro/redun/pull/217>`_ - show keyword arguments
* `216 <https://github.com/insitro/redun/pull/216>`_ - Fix isort line length
* `215 <https://github.com/insitro/redun/pull/215>`_ - DE-2623 Dont use recursive for getting execution jobs
* `213 <https://github.com/insitro/redun/pull/213>`_ - fix path term parsing
* `212 <https://github.com/insitro/redun/pull/212>`_ - fix: redun server ECS service name in merge spec
* `208 <https://github.com/insitro/redun/pull/208>`_ - Scope redun_server DB sessions at the request level
* `210 <https://github.com/insitro/redun/pull/210>`_ - Cleanup logging of migrations
* `211 <https://github.com/insitro/redun/pull/211>`_ - DE-2599 Use wait_until in aws batch tests to fix flaky tests

0.4.12
======
:Date: May 07, 2021

* `206 <https://github.com/insitro/redun/pull/206>`_ - Add method to clone RedunBackendDB with connection pool sharing
* `196 <https://github.com/insitro/redun/pull/196>`_ - DE-2325 Add database versioning commands
* `201 <https://github.com/insitro/redun/pull/201>`_ - Add quick script to generate release notes

0.4.11
======
:Date: April 22th, 2021

* `198 <https://github.com/insitro/redun/pull/198>`_ - Add support for configuration only task args
* `197 <https://github.com/insitro/redun/pull/197>`_ - [DE-2428] Fix typed list check
* `192 <https://github.com/insitro/redun/pull/192>`_ - DE-2434 Add more common tasks to functools
* `194 <https://github.com/insitro/redun/pull/194>`_ - decouple scheduler from oneshot
* `186 <https://github.com/insitro/redun/pull/186>`_ - Dockerize redun server, update directory layout and utils, add specs for prod deployment
* `190 <https://github.com/insitro/redun/pull/190>`_ - DE-2464 Add postmortem debugging

0.4.10
======
:Date: April 12th, 2021

* `188 <https://github.com/insitro/redun/pull/188>`_ - Don't let docker change terminal to raw mode
* `187 <https://github.com/insitro/redun/pull/187>`_ - Tasks should allow novel kwargs
* `180 <https://github.com/insitro/redun/pull/180>`_ - Use amazonlinux default pythons
* `185 <https://github.com/insitro/redun/pull/185>`_ - Support job timeouts on batch
* `182 <https://github.com/insitro/redun/pull/182>`_ - Lazy operators for redun Expressions

0.4.9
=====
:Date: March 23rd, 2021

* `183 <https://github.com/insitro/redun/pull/183>`_ - add py.typed
* `177 <https://github.com/insitro/redun/pull/177>`_ - Support list args from cli
* `178 <https://github.com/insitro/redun/pull/178>`_ - Fix settrace monkeypatch to restore debugging ability
* `179 <https://github.com/insitro/redun/pull/179>`_ - DE-2370 Give array jobs a unique uuid
* `181 <https://github.com/insitro/redun/pull/181>`_ - sqlalchemy 1.4.0 no longer allows postgres:// gotta be postgresql://
* `176 <https://github.com/insitro/redun/pull/176>`_ - Improve pickle preview for constructor and __new__
* `173 <https://github.com/insitro/redun/pull/173>`_ - Allow pycharm's debugger to work with redun
* `175 <https://github.com/insitro/redun/pull/175>`_ - Set choices on parser for enum args
* `174 <https://github.com/insitro/redun/pull/174>`_ - Allow use of id prefixes with push/pull commands
* `171 <https://github.com/insitro/redun/pull/171>`_ - Make S3 repositories work
* `172 <https://github.com/insitro/redun/pull/172>`_ - Match python 3.7 and 3.8 micro versions to match codebuild image


0.4.8
=====
:Date: March 10th, 2021

* `111 <https://github.com/insitro/redun/pull/111>`_ - Add concept of remote repos
* `169 <https://github.com/insitro/redun/pull/169>`_ - Remove invalid positional arg in get_or_create_job_definition call
* `147 <https://github.com/insitro/redun/pull/147>`_ - Dir should have File as subvalues for better dataflow recording
* `165 <https://github.com/insitro/redun/pull/165>`_ - Fix lack of caching for catch expressions
* `164 <https://github.com/insitro/redun/pull/164>`_ - Fix PartialTask's options() and partial() calls so that they interact correctly
* `163 <https://github.com/insitro/redun/pull/163>`_ - Imports executors in the __init__
* `155 <https://github.com/insitro/redun/pull/155>`_ - Use config_dir with redun_server

0.4.7
=====
:Date: February 24th, 2021

**WARNING:** This version contains a bug in the `get_or_create_job_defintion` call in `batch_submit`. Do not use this version.

* `156 <https://github.com/insitro/redun/pull/156>`_, `157 <https://github.com/insitro/redun/pull/157>`_, `158 <https://github.com/insitro/redun/pull/158>`_, `160 <https://github.com/insitro/redun/pull/160>`_ - Automatic publishing of packages and docs
* `153 <https://github.com/insitro/redun/pull/153>`_ - Use existing job def
* `116 <https://github.com/insitro/redun/pull/116>`_ - Display dataflow
* `154 <https://github.com/insitro/redun/pull/154>`_ - Fix data provenance recording for seq scheduler task
* `152 <https://github.com/insitro/redun/pull/152>`_ - Fix pickling expression upstreams
* `136 <https://github.com/insitro/redun/pull/136>`_ - Add redux to redun_server
* `151 <https://github.com/insitro/redun/pull/151>`_ - Record stderr from scripts on batch
* `149 <https://github.com/insitro/redun/pull/149>`_ - Add support for generating DB URI from AWS secret
* `150 <https://github.com/insitro/redun/pull/150>`_ - Document max value size
* `146 <https://github.com/insitro/redun/pull/146>`_ - Cryptic error for large falues
* `148 <https://github.com/insitro/redun/pull/148>`_ - Simplify Scheduler.run() to take expressions
* `145 <https://github.com/insitro/redun/pull/145>`_ - Add nout task option for tuples
* `144 <https://github.com/insitro/redun/pull/144>`_ - Increase sqlalchemy requirement to 1.3.17
* `143 <https://github.com/insitro/redun/pull/143>`_ - Package on submit not start

0.4.6
=====
:Date: February 3rd, 2021

* `141 <https://github.com/insitro/redun/pull/141>`_ - Only gather inflight jobs on batch on first submission

0.4.5
=====
:Date: January 28th, 2021

* `139 <https://github.com/insitro/redun/pull/139>`_ - Propagate batch script errors
* `137 <https://github.com/insitro/redun/pull/137>`_ - Override CannotInspectContainerError batch errors
* `138 <https://github.com/insitro/redun/pull/138>`_ - Fix pickle preview for classes where the module can't be found
* `133 <https://github.com/insitro/redun/pull/133>`_ - Small fixes from demo talk
* `132 <https://github.com/insitro/redun/pull/132>`_ - Small improvements to File.copy_to and self-stagin

0.4.4
=====
:Date: January 15th, 2021

* `131 <https://github.com/insitro/redun/pull/131>`_ - Fix catch dataflow
* `134 <https://github.com/insitro/redun/pull/134>`_ - Add notebook example of redun scheduler evaluation
* `128 <https://github.com/insitro/redun/pull/128>`_ - Make redun compatible with sqlalchemy-1.4.0b1
* `129 <https://github.com/insitro/redun/pull/129>`_ - Add pickle_preview for unknown classes
* `130 <https://github.com/insitro/redun/pull/130>`_ - Fix catch dataflow
* `127 <https://github.com/insitro/redun/pull/127>`_ - Add FAQ page to docs
* `126 <https://github.com/insitro/redun/pull/126>`_ - Require sorted imports

0.4.3
======
:Date: January 5th, 2021

* `122 <https://github.com/insitro/redun/pull/122>`_ - Stronger type checking for task calls
* `101 <https://github.com/insitro/redun/pull/101>`_ - Record CallNodes when an exception is raised
* `86 <https://github.com/insitro/redun/pull/86>`_ - Scheduler tasks

0.4.2
======
:Date: January 4th, 2021

* `121 <https://github.com/insitro/redun/pull/121>`_ - Array job reuniting fix

0.4.1
======
:Date: December 23rd, 2020

* `119 <https://github.com/insitro/redun/pull/119>`_ - Bugfix to correctly restart job array monitor thread

0.4.0
======
:Date: December 15th, 2020

* `83 <https://github.com/insitro/redun/pull/83>`_ - Detect and submit job arrays to AWS batch
* `114 <https://github.com/insitro/redun/pull/114>`_ - Adds job definition option to run container in privileged mode

0.3.12
======
:Date: December 10th, 2020

* `76 <https://github.com/insitro/redun/pull/76>`_ - Improve querying of logs

0.3.11
======
:Date: December 8th, 2020

* `109 <https://github.com/insitro/redun/pull/109>`_ - Permalink update in README
* `108 <https://github.com/insitro/redun/pull/108>`_ - Automated release

0.3.10
======
:Date: December 3rd, 2020

* `104 <https://github.com/insitro/redun/pull/104>`_ - use ECR for postgres image
* `95 <https://github.com/insitro/redun/pull/95>`_ - Hard fail on script errors
* `100 <https://github.com/insitro/redun/pull/100>`_ - Show more information in logs and traceback
* `102 <https://github.com/insitro/redun/pull/102>`_ - Fix check-valid=shallow to use the original call node
* `98 <https://github.com/insitro/redun/pull/98>`_ - Skip license check when building conda packages
* `105 <https://github.com/insitro/redun/pull/105>`_ - Typecheck map_nested_value
* `103 <https://github.com/insitro/redun/pull/103>`_ - Fix script reactivity to inputs and outputs
* `106 <https://github.com/insitro/redun/pull/106>`_ - Small clean up of batch logs

0.3.9
=====
:Date: November 25th, 2020

* `96 <https://github.com/insitro/redun/pull/96>`_ - Default to interactive debugging
* `81 <https://github.com/insitro/redun/pull/81>`_ - Allow REDUN_CONFIG environment variable to specify config directory
* `92 <https://github.com/insitro/redun/pull/92>`_ - DE-1922 tolerate missing logs for failed jobs

0.3.8
=====
:Date: November 18th, 2020

* `89 <https://github.com/insitro/redun/pull/89>`_ - Respect no-cache for job reuniting.
* `88 <https://github.com/insitro/redun/pull/88>`_ - Assume batch output after completion is valid.
* `87 <https://github.com/insitro/redun/pull/87>`_ - Fix filesystem caching and Dir hashing caching.
* `85 <https://github.com/insitro/redun/pull/85>`_ - Add step to publish pypi package in publish script.
* `84 <https://github.com/insitro/redun/pull/84>`_ - Fix package name in dependencies notes in README.

0.3.7
=====
:Date: November 12th, 2020

* `80 <https://github.com/insitro/redun/pull/80>`_ - redun import paths should take precedence over system imports.
* `79 <https://github.com/insitro/redun/pull/79>`_ - fix default arg parsing and prefix args.

0.3.6
=====
:Date: November 10th, 2020

* `73 <https://github.com/insitro/redun/pull/73>`_ - Allow users to customize `setup_scheduler()`.

0.3.5
=====
:Date: November 10, 2020

* `77 <https://github.com/insitro/redun/pull/77>`_ - Check version of redun cli in docker container.

0.3.4
=====
:Date: October 29th, 2020

* `72 <https://github.com/insitro/redun/pull/72>`_ - Use current working directory when importing a module.
* `64 <https://github.com/insitro/redun/pull/64>`_ - Some optimizations for AWS Batch large fanout.  

0.3.3
=====
:Date: October 28th, 2020

* `#71 <https://github.com/insitro/redun/pull/71>`_ - Don't fetch batch logs when debug=True

0.3.2
=====
:Date: October 27th, 2020

* `#66 <https://github.com/insitro/redun/pull/66>`_ - Fix import_script to properly support module-style

0.3.1
=====

* Fix bug with using s3fs >= 0.5

0.3
=====
:Date: October 20th, 2020

* Improve display of errors and logs for AWS Batch jobs.

0.2.5
=====
:Date: October 14th, 2020

* `#57 <https://github.com/insitro/redun/pull/57>`_ - Improve redun traceback for failed jobs.
* `#56 <https://github.com/insitro/redun/pull/56>`_ - Fix local shell error propogation.
* `#54 <https://github.com/insitro/redun/pull/54>`_ - Add documentation on required dependencies.

0.2.4
=====
:Date: October 6, 2020

* Encourage defining task namespaces by raising a warning. The warning can be ignored using a [configuration option](config.html#ignore-warnings).


0.2.3
=====
:Date: September 25, 2020

* Fixes FileNotFoundError occuring when using AWS Batch tasks, by avoiding the s3fs cache.


0.2.2
=====
:Date: August 27, 2020

* Require database credentials to be specified by environment variables


0.2.1
=====

:Date: August 9, 2020

 * Fix duplicate upstream bug.


0.2.0
=====

:Date: August 7, 2020

 * Add support for Python 3.8


0.1.1
=====

:Date: July 29, 2020

 * Drop dependency on bcode as it has no conda package and the repo appears abandoned.


0.1
===

 * Initial release.
