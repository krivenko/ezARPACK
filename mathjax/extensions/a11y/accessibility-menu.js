/*
 *  /MathJax/extensions/a11y/accessibility-menu.js
 *
 *  Copyright (c) 2009-2018 The MathJax Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

!function(g,h){var p,k,j=g.config.menuSettings,m=Function.prototype.bind?function(i,a){return i.bind(a)}:function(i,a){return function(){i.apply(a,arguments)}},c=Object.keys||function(i){var a=[];for(var l in i){i.hasOwnProperty(l)&&a.push(l)}return a},d=MathJax.Ajax.config.path;d.a11y||(d.a11y=g.config.root+"/extensions/a11y");var f=h["accessibility-menu"]={version:"1.6.0",prefix:"",defaults:{},modules:[],MakeOption:function(a){return f.prefix+a},GetOption:function(a){return j[f.MakeOption(a)]},AddDefaults:function(){for(var o,l=c(f.defaults),q=0;o=l[q];q++){var a=f.MakeOption(o);void 0===j[a]&&(j[a]=f.defaults[o])}},AddMenu:function(){for(var s,r=Array(this.modules.length),v=0;s=this.modules[v];v++){r[v]=s.placeHolder}var q,l,u=k.FindId("Accessibility");u?(r.unshift(p.RULE()),u.submenu.items.push.apply(u.submenu.items,r)):((q=(k.FindId("Settings","Renderer")||{}).submenu)&&(r.unshift(p.RULE()),r.unshift(q.items.pop()),r.unshift(q.items.pop())),r.unshift("Accessibility"),u=p.SUBMENU.apply(p.SUBMENU,r),(l=k.IndexOfId("Locale"))?k.items.splice(l,0,u):k.items.push(p.RULE(),u))},Register:function(a){f.defaults[a.option]=!1,f.modules.push(a)},Startup:function(){p=MathJax.Menu.ITEM,k=MathJax.Menu.menu;for(var i,a=0;i=this.modules[a];a++){i.CreateMenu()}this.AddMenu()},LoadExtensions:function(){for(var i,a=[],l=0;i=this.modules[l];l++){j[i.option]&&a.push(i.module)}return a.length?g.Startup.loadArray(a):null}},b=MathJax.Extension.ModuleLoader=MathJax.Object.Subclass({option:"",name:["",""],module:"",placeHolder:null,submenu:!1,extension:null,Init:function(r,q,s,o,l){this.option=r,this.name=[q.replace(/ /g,""),q],this.module=s,this.extension=o,this.submenu=l||!1},CreateMenu:function(){var a=m(this.Load,this);this.submenu?this.placeHolder=p.SUBMENU(this.name,p.CHECKBOX(["Activate","Activate"],f.MakeOption(this.option),{action:a}),p.RULE(),p.COMMAND(["OptionsWhenActive","(Options when Active)"],null,{disabled:!0})):this.placeHolder=p.CHECKBOX(this.name,f.MakeOption(this.option),{action:a})},Load:function(){g.Queue(["Require",MathJax.Ajax,this.module,["Enable",this]])},Enable:function(i){var a=MathJax.Extension[this.extension];a&&(a.Enable(!0,!0),MathJax.Menu.saveCookie())}});f.Register(b("collapsible","Collapsible Math","[a11y]/collapsible.js","collapsible")),f.Register(b("autocollapse","Auto Collapse","[a11y]/auto-collapse.js","auto-collapse")),f.Register(b("explorer","Explorer","[a11y]/explorer.js","explorer",!0)),f.AddDefaults(),g.Register.StartupHook("End Extensions",function(){g.Register.StartupHook("MathMenu Ready",function(){f.Startup(),g.Startup.signal.Post("Accessibility Menu Ready")},5)},5),MathJax.Hub.Register.StartupHook("End Cookie",function(){MathJax.Callback.Queue(["LoadExtensions",f],["loadComplete",MathJax.Ajax,"[a11y]/accessibility-menu.js"])})}(MathJax.Hub,MathJax.Extension);
