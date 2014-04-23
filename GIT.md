git
===


Install
-------

plain user:

 * in terminal execute: 
 
    `cd basedirectory`
 
    `git clone https://github.com/ChristophKirst/StemCellTracker`


developer / programmer:

  * create account at [http://github.com](http://github.com)

  * got to [https://github.com/ChristophKirst/StemCellTracker](https://github.com/ChristophKirst/StemCellTracker) and press fork button 

  * in terminal execute:
	
	`cd basedirectory`

	`git clone https://github.com/username/StemCellTracker.git`
	
  * configure remotes (named upstream)
        
	`cd StemCellTracker`

	`git remote add upstream https://github.com/ChristophKirst/StemCellTracker`

	`git fetch upstream`


Backup
------

to backup your version in case you followed the developer / proogrammer route:

  * in termianl in the StemCellTracker directory execute:

      `git add -A`

      `git commit -m 'some description of what you did'`

      `git push`


Update
------    

plain user:

  * in terminal in the StemCellTracker directory execute
     
      `git pull`


programmer: 

in case you want to update your code from the upstream repository

  * in terminal execute:
 
      `git fetch upstream`
      
      `git merge upstream/master`

  * if mergin fails, some files will be highlighted with <<<<<< >>>>>> entries, fix this manually

  * if you dont care about your own changes and simply want the plain new version:

      `git reset --hard upstream/master`

  * to force it to your fork on github use
       
	  `git push origin/master --force` 


Submitting
----------

in case you have something to contribute to the code:
 
  * follow the steps in the Backup section first

  * got to [http://github.com](http://github.com) and click pull request 


Refs
----

stem cell tracker home:

  * [https://github.com/ChristophKirst/StemCellTracker](https://github.com/ChristophKirst/StemCellTracker)

A good source to get questions answered: 

  * [https://help.github.com/](https://help.github.com)
  * [http://git-scm.com/documentation](http://git-scm.com/documentation)

git home:

  * [http://git-scm.com](http://git-scm.com/)

