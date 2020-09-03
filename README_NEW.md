# ANUGA Community TechLauncher Project

## Overview

The existing ANUGA Hydro software has been used to model floods for more than a decade, with excellent results. But for smaller floods, it has a tendency to over estimate the flood area, due to being unable to model underground drainage systems.

This project will extend the ANUGA Hydro software, which is capable of hydrodynamic modelling, by coupling with the US EPA's Storm Water Management Model (SWMM), thus adding to it the ability to model the effects of underground drainage systems. 

## Team

|  Name          | UID    | Principal Role | Secondary Role | Other |
|:--------------:|:------:|:--------------:|:--------------:|:-----:|
| Zhixian (Rachel) Wu | u5807060 |  |  | Spokesperson, Team Leader |
| Xingnan Pan | u6744662 |  |  | Deputy Spokesperson, Project Manager |
| Chen Huang | u6735118 |  |  |  |
| Lixin Hou | u6456457 |  |  |  |
| Mingda Zheng | u6686733 |  |  |  |
| Yijie Liu | u6890141 |  |  |  |
| Zijun Zhang | u6904534 |  |  |  |

## Stakeholders
* **The sponsors:**
   * Professor Stephen Roberts, ANU
   * Dr Ole Nielsen, Geoscience Australia
* **The user representatives (flood modellers):**
   * Rudy Van Drie, Balance Research and Development
   * Dr Petar Milevski, Civil Engineer Urban Drainage, Wollongong City Council

## Documentation

### Sprint Stories

> [Trello](https://trello.com/b/Z45C7crP/agile-sprint-board)

### Communication

> [Slack](https://anu-flood-modelling.slack.com)

### Meeting Minutes

Sprint 1 (start of semester - 19/08/2020)

> [2020-08-04 Team Meeting](https://docs.google.com/document/d/1SW3PUsRs-9bc1CYlVkW6fHQLiOQ0cm0w_jzSKu37Gpw/edit?usp=sharing)
> [2020-08-06 Client Meeting](https://docs.google.com/document/d/1J_kqxAhOHSAh3xWV8enVu0XkZSba1jQchf01azwkgvg/edit?usp=sharing)

### Decisions

> [Log for Small Decisions](https://docs.google.com/spreadsheets/d/1uPZlRMNaRBlZnUdfNPVQ4e_S48npiRRkqP9GHJUyXS4/edit?usp=sharing)
> [Documents for Large Decisions](https://docs.google.com/spreadsheets/d/1uPZlRMNaRBlZnUdfNPVQ4e_S48npiRRkqP9GHJUyXS4/edit?usp=sharing)

## Timeline

We are doing two-week sprints, with client meetings to close each sprint on Wednesday 5:00PM Canberra time, and team meetings for sprint retrospectives and sprint planning on Wednesday 7:00PM Canberra time.

The first sprint will be a bit longer, so that the rest of the sprints will end just before the Week 6 and Week 10 audits. This means the first sprint will end Wednesday of Week 4.

* **2020-s2-timeline**
<img src="https://drive.google.com/file/d/1fBOS3L8SHISkuszNiWiaJWVJX4TLfvgx/view?usp=sharing" alt="20-s2-TL" align=center />

* **2021-s1-timeline (TBC)**

## Risks

|Risk ID|Risk points|Mitigation measures|
|:-----:|:---------:|:-----------------:|
|1|The time difference might be a cooperation barrier as the team consists of offshore and native members|Most offshore members are living in China, which merely has 2 hours lag with the Australian standard time. Therefore, the team or client meeting can be set at afternoon to mitigate the impact.|
|2|It may take longer time for team members to learn the complex models in the project, as most members are new to machine learning area or lack of experience.|The team is able to split each stage task into some small tasks which will be allocated to a small group of members (i.e. 1~3 members). Each small tasks will be learnt and conducted simultaneously|
|3|Team members may have some emergencies during the project, such as sick, exam, which may interrupt the project progress.|Each task will be allocated with task supervisor(s) and a shadow group. The supervisor(s) will conduct the task and the shadow  group should monitor the progress and checkout the performance. If the supervisor(s) occurs some emergency leading to lack of members, the shadow team is able to embed into the progress immediately.|
|4|The system and equipment requirements may cause some difficulties to the team, as the project is required to design in Ubuntu 20.04, but some features can only run in Ubuntu 18.04. In addtion, the models require high capacity RAM to get a good performance.|Members can use virtual machine or install dual systems to match the development circumstance. And some complex issues can be tested in lab machines by members in Canberra.| 

## Technical Constraints

* The project should be developed in Github
   * Each member is able to test in a branch
   * Using pull request to get the task review from others
   * Merging codes into the `main` branch

* The project is mainly conducted on Ubuntu 20.04

* Setup Continuous Integration (CI) tools to test on three platfroms (Windows, MacOS and Ubuntu) automatically

* The project should follow the Agile process 
   * Two week length sprint
   * Create backlog ([Trello](https://trello.com/b/Z45C7crP/agile-sprint-board)) to track issues

* Software standards
   * The Python code should follow the [PEP8](https://www.python.org/dev/peps/pep-0008/) standard apart from agreed exceptions
   * All codes should have corresponding unit tests apart from the most trivial
   * Model behavior should be tested end to end
   * All tests should be included in test suites that can be run either manually or by the CI server.
