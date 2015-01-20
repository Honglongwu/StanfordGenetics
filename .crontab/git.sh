#crontab -e
# */30 * * * * bash /srv/gsfs0/projects/steinmetz/hansun/.crontab/git.sh
cd /srv/gsfs0/projects/steinmetz/hansun
git add *.py
git add *.sh
git add *.R
git commit -a -m 'hanice'
git push origin master
