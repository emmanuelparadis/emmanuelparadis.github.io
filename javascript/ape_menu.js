function openNav() {
document.getElementById("mySidebar").style.width="250px";
document.getElementById("main").style.marginLeft="250px";
}

function closeNav() {
document.getElementById("mySidebar").style.width="0";
document.getElementById("main").style.marginLeft="0";
}

changeDisplayState = function (id, id2, str1, str2) {
     var d = document.getElementById(id2),
          e = document.getElementById(id);
     if (e.style.display === 'none' || e.style.display === '') {
          e.style.display = 'block';
          d.innerHTML = str2;
     }
     else {
          e.style.display = 'none';
          d.innerHTML = str1;
     }
};

