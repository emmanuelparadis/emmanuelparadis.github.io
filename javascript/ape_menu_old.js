var timeout = 250;
var closetimer = 0;
var ddmenuitem = 0;

// deroule un menu
function mopen(id)
{
    mcancelclosetime();

    // referme le menu precedent
    if (ddmenuitem) ddmenuitem.style.visibility = 'hidden';

    // deroule le nouveau menu
    ddmenuitem = document.getElementById(id);
    ddmenuitem.style.visibility = 'visible';

}

// ferme un menu
function mclose()
{
    if (ddmenuitem) ddmenuitem.style.visibility = 'hidden';
}

function mclosetime()
{
    closetimer = window.setTimeout(mclose, timeout);
}

function mcancelclosetime()
{
    if (closetimer)
	{
	    window.clearTimeout(closetimer);
	    closetimer = null;
	}
}

// referme un menu quand la souris n'est plus dessus
document.onclick = mclose;

