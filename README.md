# ANUGA Community TechLauncher Project

## Table of Contents

1. [Overview](#overview)
2. [Team](#team)
3. [Stakeholders](#stakeholders)
4. [Documentation](#documentation)
5. [Timeline](#timeline)
6. [Risks](#risks)
7. [Tools and Client Requirements](#tools-and-client-requirements)
8. [Technical Constraints](#technical-constraints)

## Overview

The existing ANUGA Hydro software has been used to model floods for more than a decade, with excellent results. But for smaller floods, it has a tendency to over estimate the flood area, due to being unable to model underground drainage systems.

This project will extend the ANUGA Hydro software, which is capable of hydrodynamic modelling, by coupling with the US EPA's Storm Water Management Model (SWMM), thus adding to it the ability to model the effects of underground drainage systems. 

## Team

|  Name          | UID    | Principal Role | Secondary Role |
|:--------------:|:------:|:---------------|:---------------|
| Zhixian (Rachel) Wu | u5807060 | SWMM developer<br>Team leader<br>Minutes taker<br>Landing page maintainer | Documentation (grammar) reviewer<br>Code reviewer (ANUGA)<br>Client manager |
| Xingnan Pan | u6744662 | ANUGA developer<br>Project Manager | Test developer<br>Code reviewer (SWMM) |
| Chen Huang | u6735118 | SWMM developer<br>Documention author<br>Minutes taker | Code reviewer (ANUGA)<br>Quality manager<br>Fallback project manager |
| Lixin Hou | u6456457 | ANUGA developer<br>Documantation author | Client Manager |
| Mingda Zheng | u6686733 | ANUGA developer<br>Quality manager | Code reviewer (SWMM)<br>Test developer|
| Yijie Liu | u6890141 | ANUGA developer<br>Quality manager<br>Documentation author | Test developer<br> |
| Zijun Zhang | u6904534 | SWMM developer<br>Documantation author | Code reviewer (ANUGA) |

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

> [2020-08-06 Team Meeting](https://docs.google.com/document/d/1J_kqxAhOHSAh3xWV8enVu0XkZSba1jQchf01azwkgvg/edit?usp=sharing)

> [2020-08-11 Team Meeting](https://docs.google.com/document/d/1Y9Rpm179KeMyxmN4nzFfLMplG6DxY44MxcMesnc1g28/edit?usp=sharing)

### Decisions

> [Log for Small Decisions](https://docs.google.com/spreadsheets/d/1uPZlRMNaRBlZnUdfNPVQ4e_S48npiRRkqP9GHJUyXS4/edit?usp=sharing)

> [Template for Large Decisions](https://docs.google.com/document/d/11qM3a2_Abr2oGtYLgIPA4QjgSELj9RFD4IboVFuBqEg/edit?usp=sharing)

> [2020-08-05 Continuous Integration Tool Selection](https://docs.google.com/document/d/1xt46NBabq5xelkVywf4NLt9Su33GicAFldKurt747fs/edit?usp=sharing)

## Timeline

We are doing two-week sprints, with client meetings to close each sprint on Wednesday 5:00PM Canberra time, and team meetings for sprint retrospectives and sprint planning on Wednesday 7:00PM Canberra time.

The first sprint will be a bit longer, so that the rest of the sprints will end just before the Week 6 and Week 10 audits. This means the first sprint will end Wednesday of Week 4.

* **2020-s2-timeline**
<img src="https://github.com/rachelwu21/anuga_core/blob/master/20-s2%20timeline.jpg" alt="20-s2-TL" align=center />

* **2021-s1-timeline (TBC)**

## Risks

|Risk ID|Risk points|Mitigation measures|
|:-----:|:---------:|:-----------------:|
|1|The time difference might be a cooperation barrier as the team consists of overseas and native members|Most overseas members are living in China, which merely has 2 hours lag with the Australian Eastern Standard Time. The remaining member lives in Perth, which is in the same timezone as China. Therefore, the team or client meeting can be set at afternoon to mitigate the impact.|
|2|It may take a long time for team members to learn the complex models in the project.|The team is able to split each stage task into some small tasks which will be allocated to a small group of members (i.e. 1~3 members). Each small tasks will be learnt and conducted simultaneously. The members can then write short sprint reports to document and report their progress and discoveries to the other members.|
|3|Team members may have some emergencies during the project, such as sick, exam, which may interrupt the project progress.|We will never have any task that is only performed by one team member. Either the task will be performed by a small group, or if it is too small one team member will be assigned as the secondary person responsible for reviewing the code and taking over if the member principally responsible has an emergency situation.|

## Tools and Client Requirements

* The project should be developed in Github
   * Each member is able to test in a branch
   * Using pull request to get the task review from others
   * Only tested and review code should be merged into the `main` branch
* The project is mainly developed on Ubuntu 20.04
   * This means that team members will need to install a virtual machine or dual boot. All members have already done so.
* Setup Continuous Integration (CI) tools to test on three platforms (Windows, MacOS and Ubuntu) automatically.
   * This was a Sprint 1 task for two members of the team. They have already set up Appveyor and TravisCI to handle this.
* Software standards
   * The Python code should follow the [PEP8](https://www.python.org/dev/peps/pep-0008/) standard apart from agreed exceptions.
   * All code, apart from the most trivial, should have corresponding unit tests.
   * Model behaviour should be tested end to end with real data.
   * Tests should be integrated with a CI server.
* The standard official version of SWMM from the US EPA website is only available for Windows, so we will use another open-source project called PySWMM by Open Water Analytics.

## Technical Constraints

The end modelling software must be a coupling between ANUGA and SWMM. There are no other open-source options for this type of software. And even if there were, the team was commissioned by the clients to improve the existing ANUGA Hydro software in a specific way. 

## Build Status
[![Build Status](https://travis-ci.com/20-S2-2-C-Flood-Modelling/anuga_core.svg?branch=master)](https://travis-ci.com/20-S2-2-C-Flood-Modelling/anuga_core)

